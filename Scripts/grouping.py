from typing import Any, Dict, List, Optional, Sequence
import scipy.stats as stats # type: ignore
from grouping_model import CSInfo,  Var, CSLocus, Locus, GroupingOptions, Grouping, PeakLocus, LDMode
from load_tabix import TabixResource
from data_access.db import Variant,CS, CSAccess, CSVariant, LDAccess
from data_access.csfactory import csfactory


def ld_threshold(ld_thresh:float, mode: LDMode,pval:float)->float:
    if mode == LDMode.DYNAMIC:
        return float(min(ld_thresh/stats.chi2.isf(pval,df=1),1.0))
    else:
        return ld_thresh

def credible_grouping(credsets:List[CS],summstat_resource: TabixResource,ld_api: LDAccess, options: GroupingOptions) -> Dict[Variant,List[Var]]:
    """Return ld partners for credible sets. Returned as a dictionary of lead var -> ld partner keys & values
    """
    output = {}
    #create list of credible set lead variants
    cs_lead_vars = [a.lead for a in credsets]
    #load pvals from summstat
    cs_pvals = {}
    for c in cs_lead_vars:
        r_c = c.chrom
        r_s = c.pos
        r_e = c.pos+1
        d = summstat_resource.load_region(r_c,r_s,r_e,[options.column_names.pval])
        try:
            cs_pvals[c] = float(d[c][0])
        except:
            raise Exception(f"Fatal Error: Credible set lead variant {c} was not found in the summary statistic resource. Can not continue.")
    csleads_by_pval = sorted([(key,value)for key,value in cs_pvals.items()],key=lambda x:x[1])
    variants_that_can_not_be_included:set[Variant] = set()
    for lead_var,lead_pval in csleads_by_pval:
        # get LD partners
        #dynamic/static r2 thresh decision
        ld_thresh = ld_threshold(options.r2_threshold,options.ld_mode,lead_pval)
        print(lead_var,lead_pval,ld_thresh)
        ld_partners = ld_api.get_range(lead_var,options.range,ld_thresh)
        ld_partners = [a for a in ld_partners if a.variant1 == lead_var]
        # get locus in summary statistic
        summstat_data = {a:{"pval":float(b[0]),"beta":float(b[1])} for a,b in summstat_resource.load_region(lead_var.chrom,lead_var.pos-options.range,lead_var.pos+options.range,[
            options.column_names.pval,options.column_names.beta
        ]).items()}
        #filter by p-value
        summstat_data = {a:b for a,b in summstat_data.items() if b["pval"]<= options.p2_threshold}
        #filter ld partners by summstat_data.
        #   Easiest option is using in dict keys, if not the most performant. 
        #   However, I don't know the performance of these different options -> just do it one way and worry about it later.
        ld_partners = [a for a in ld_partners if a.variant2 in summstat_data]
        #from LD partners, exclude all current cs members, all cs leads, and any other variants.
        current_cs = [a for a in credsets if a.lead == lead_var][0]
        filter_list = [current_cs.lead]
        filter_list.extend([a.variant for a in current_cs.variants])
        filter_set = set(filter_list)
        ld_partners = [a for a in ld_partners if ((a.variant2 not in filter_set)and (a.variant2 not in variants_that_can_not_be_included))]

        
        #build the ld partner data
        ld_partner_data = [
            Var(a.variant2,
                summstat_data[a.variant2]["pval"],
                summstat_data[a.variant2]["beta"],
                a.r2)
            for a in ld_partners
        ]
        #variants_that_can_not_be_included.add(lead_var)
        if not options.overlap:
            variants_that_can_not_be_included = variants_that_can_not_be_included.union([lead_var])
            variants_that_can_not_be_included = variants_that_can_not_be_included.union([a.variant for a in current_cs.variants])
            variants_that_can_not_be_included = variants_that_can_not_be_included.union([a.id for a in ld_partner_data])
        output[lead_var] = ld_partner_data
        # add to ld partners, filter out those variants that can not be ld partners anymore (ld partners, cs variants)
        #basically iterate over the leads in ascending pval order, gathering the LD partners.
    #order by ascending p-value
    #for each cs, by their lead variant
    #add ld variants. Include also other CS variants if they are included, and mark them like they were in the original implementation. Filter out variants in the credset.
    #keep a list of variants that should not be included in the future ld variants if overlap is false
    return output

def in_range(v1:Var, v2:Var,pos_range:int)->bool:
    """Return True if |v1-v2| <= range, otherwise False
    """
    if (v1.id.chrom == v2.id.chrom) and (v1.id.pos <= v2.id.pos + pos_range ) and (v1.id.pos >= v2.id.pos - pos_range):
            return True
    return False

def simple_grouping(summstat_resource:TabixResource, options:GroupingOptions) -> List[PeakLocus]:
    """Form autoreporting groups with simple, range-based grouping
    """
    # sift through the data, piling it into 1) pile of possible leads, 2) pile of range-based partners.
    # order by p-value the first, probabaly best to also divide into chromosomes at this point
    # start with the highest p-value, group all those variants in either pile that are close enough to it
    # then remove from consideration, probably best to keep a set of already grouped, so as not to modify data.
    output = []
    p1_pile:List[Var] = []
    p2_pile:List[Var] = []
    cpra = [
        options.column_names.c,
        options.column_names.p,
        options.column_names.r,
        options.column_names.a,
        options.column_names.pval,
        options.column_names.beta
    ]
    ## load p1 and p2 filtered variants from summary statistic
    hd = summstat_resource.header
    hdi = {a:i for i,a in enumerate(hd)}
    for l in summstat_resource.fileobject.fetch():
        cols = l.split("\t")
        pval = float(cols[hdi[cpra[4]]])
        if pval < options.p2_threshold:
            beta = float(cols[hdi[cpra[5]]])
            var = Var(Variant(
                cols[hdi[cpra[0]]],
                int(cols[hdi[cpra[1]]]),
                cols[hdi[cpra[2]]],
                cols[hdi[cpra[3]]]),
                pval,
                beta,
                None
            )
            p2_pile.append(var)
            if pval < options.p1_threshold:
                p1_pile.append(var)
    p1_pile = sorted(p1_pile,key=lambda x: (x.pval,x.id.chrom,x.id.pos),reverse=True)
    while p1_pile:
        lead_var = p1_pile.pop()
        p2_rangevars = [a for a in p2_pile if (in_range(lead_var,a,options.range) and a != lead_var)]
        remove_set = set(lead_var)
        remove_set = remove_set.union(p2_rangevars)
        #filter piles
        p1_pile = [a for a in p1_pile if a not in remove_set]
        if not options.overlap:
            p2_pile = [a for a in p2_pile if a not in remove_set]
        else:
            p2_pile = [a for a in p2_pile if a != lead_var]
        output.append(PeakLocus(lead_var,p2_rangevars,Grouping.RANGE))
    #oh shit this need to be completed
    return output

def ld_grouping(summstat_resource: TabixResource,ld_api:LDAccess,options:GroupingOptions)->List[PeakLocus]:
    """Form autoreporting groups with LD grouping
    """
    # sift the dta to create the candidate variants in either group. Probably should not have duplicates in the groups.
    # Then, going from most significant to least significant, start grouping
    # here, one has to join the ld output to the variants.
    # Then, keep up a list of variants that can not be added to groups since they have been eliminated.
    # in that case, actually should have all p1 vars in p2 pile.
    ## Init variables
    output = []
    p1_pile:List[Var] = []
    p2_pile:List[Var] = []
    cpra = [
        options.column_names.c,
        options.column_names.p,
        options.column_names.r,
        options.column_names.a,
        options.column_names.pval,
        options.column_names.beta
    ]
    ## load p1 and p2 filtered variants from summary statistic
    hd = summstat_resource.header
    hdi = {a:i for i,a in enumerate(hd)}
    for l in summstat_resource.fileobject.fetch():
        cols = l.split("\t")
        pval = float(cols[hdi[cpra[4]]])
        if pval < options.p2_threshold:
            beta = float(cols[hdi[cpra[5]]])
            var = Var(Variant(
                cols[hdi[cpra[0]]],
                int(cols[hdi[cpra[1]]]),
                cols[hdi[cpra[2]]],
                cols[hdi[cpra[3]]]),
                pval,
                beta,
                1.0
            )
            p2_pile.append(var)
            if pval < options.p1_threshold:
                p1_pile.append(var)
    ## order p1 piel to be a stack with most significant variant at the top (i.e. last value)
    p1_pile = sorted(p1_pile,key=lambda x: x.pval,reverse=True)
    ## Grouping algorithm
    while p1_pile:
        lead_var = p1_pile.pop()
        lead_variant = lead_var.id
        #static/dynamic r2
        ld_thresh = ld_threshold(options.r2_threshold,options.ld_mode,lead_var.pval)
        ld_data = ld_api.get_range(lead_variant,options.range,ld_thresh)
        ld_data = [a for a in ld_data if (a.variant1 == lead_variant) and (a.variant2 != lead_variant)]
        ld_data_varset = set([a.variant2 for a in ld_data])
        #filter p2 pile by ld data. 
        ld_partners= [a for a in p2_pile if a.id in ld_data_varset]
        ## ld_partner_set = set(ld_partners)
        #Now, ld_partner_idscontains all of the required ld partners:
        # The definition for ld partners is 1) pval < p2_threshold, 2) in sufficient LD with the lead var
        # could also be done by loading a tabix region for each of the LD regions... oh well
        #join r2 value from ld_data to ld partners
        #first, do a index dict
        ld_index_dict = {a.variant2:i for i,a in enumerate(ld_data) }
        #after that getting the right LD value from ld data is trivial, since the LD data is filtered to be variant1 = lead variant
        ld_partners = [Var(a.id,a.pval,a.beta,ld_data[ld_index_dict[a.id]].r2) for a in ld_partners]
        #filter all p2 variants out of p1 variants. This along with the pop in the beginning ensures convergence
        p1_pile = [a for a in p1_pile if a.id not in ld_data_varset]
        #if not overlap, filter the p2_pile not to contain the already joined ld partners.
        if not options.overlap:
            p2_pile = [a for a in p2_pile if a.id not in ld_data_varset]
        ## Form locus
        output.append(PeakLocus(
            lead_var,
            ld_partners,
            Grouping.LD
        ))
    return output

def filter_gws_variants(summstat_resource: TabixResource, options: GroupingOptions)->List[Var]:
    p1_pile = []
    cpra = [
        options.column_names.c,
        options.column_names.p,
        options.column_names.r,
        options.column_names.a,
        options.column_names.pval,
        options.column_names.beta
    ]
    hd = summstat_resource.header
    hdi = {a:i for i,a in enumerate(hd)}
    for l in summstat_resource.fileobject.fetch():
        cols = l.split("\t")
        pval = float(cols[hdi[cpra[4]]])
        if pval < options.p1_threshold:
            beta = float(cols[hdi[cpra[5]]])
            var = Var(Variant(
                cols[hdi[cpra[0]]],
                int(cols[hdi[cpra[1]]]),
                cols[hdi[cpra[2]]],
                cols[hdi[cpra[3]]]),
                pval,
                beta,
                None
            )
            p1_pile.append(var)
    return p1_pile

def form_groups(summstat_resource: TabixResource, gr_opts:GroupingOptions, cs_access: CSAccess, ld_access:LDAccess) -> Sequence[Locus]:
    """Group summary statistics according to grouping options 
    """
    #load cs
    group_type = gr_opts.grouping_mode
    if group_type == Grouping.NONE:
        #filter summstat to gws variants only
        vars = filter_gws_variants(summstat_resource,gr_opts)
        loci:Sequence[Locus] = [PeakLocus(v,None,Grouping.NONE) for v in vars]
    elif group_type == Grouping.RANGE:
        #range-based grouping
        loci = simple_grouping(summstat_resource,gr_opts)
    elif group_type == Grouping.LD:
        #ld-based grouping
        loci = ld_grouping(summstat_resource,ld_access,gr_opts)
    elif group_type == Grouping.CS:
        if cs_access:
            cs = cs_access.get_cs()
        else:
            raise Exception(f"Credible set grouping mode set, but credible set file was not provided or it was not loaded correctly!")
        ### Constructing valid credible sets
        cs_variants = list(set([a.variant for credset in cs for a in credset.variants]))
        cs_pvalbetadict = {}
        for v in cs_variants:
            cs_pvalbetadict.update( {a:[float(c) for c in b] for a,b in summstat_resource.load_region(v.chrom,v.pos,v.pos+1,[gr_opts.column_names.pval,gr_opts.column_names.beta]).items()} )
        cs_groups = {}
        cs_infos = {}
        #get pval, beta, r2 for credible set variants
        for c in cs:
            #gather ld data
            c_variants = set([a.variant for a in c.variants])
            max_range = max([abs(c.lead.pos - a.variant.pos) for a in c.variants])+5000
            #susie might not have LD so calculate it here
            ld_data = ld_access.get_range(c.lead,max_range,0.0)
            ld_data = [a for a in ld_data if ( (a.variant1 == c.lead) and (a.variant2 in c_variants))  ]
            ld_dict = {a.variant2:a.r2 for a in ld_data}

            vs = [Var(
                    v.variant,
                    cs_pvalbetadict.get(v.variant,(float("nan"),float("nan")))[0],
                    cs_pvalbetadict.get(v.variant,(float("nan"),float("nan")))[1],
                    ld_dict.get(v.variant,float("nan"))
                )
                for v in c.variants]
            lead_var = Var(c.lead,cs_pvalbetadict.get(c.lead,(float("nan"),float("nan")))[0],
                    cs_pvalbetadict.get(c.lead,(float("nan"),float("nan")))[1],
                    1.0)
            cs_groups[lead_var] = vs

            #add csinfo for csinfos
            cs_info = CSInfo(
                c.region,
                c.lead,
                c.number,
                c.bayes,
                c.min_r2,
                c.size,
                c.good_cs
            )
            cs_infos[lead_var] = cs_info
        ### Get LD partners, construct loci.
        ld_partner_dict = credible_grouping(cs,summstat_resource,ld_access,gr_opts)

        loci = [CSLocus(a,cs_groups[a],cs_infos[a],ld_partner_dict[a.id]) for a in cs_groups.keys()]

    else:
        raise Exception("This grouping type not implemented!")
    
    return loci
