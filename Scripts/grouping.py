from typing import Any, Dict, List, Optional, Sequence
import scipy.stats as stats # type: ignore
from grouping_model import CSInfo,  Var, CSLocus, Locus, GroupingOptions, Grouping, PeakLocus, LDMode
from load_tabix import TabixResource
from data_access.db import Variant,CS, CSAccess, CSVariant, LDAccess
from data_access.csfactory import csfactory
from time_decorator import timefunc

def ld_threshold(ld_thresh:float, mode: LDMode,pval:float, pval_is_mlog10p:bool=False)->float:
    if mode.value == LDMode.DYNAMIC.value:
        raw_pval = max(10**(-pval), 5e-324) if pval_is_mlog10p else pval
        return float(min(ld_thresh/stats.chi2.isf(raw_pval,df=1),1.0))
    else:
        return ld_thresh

def _passes_threshold(value:float, threshold:float, pval_is_mlog10p:bool)->bool:
    """Check if value passes significance threshold.
    For raw pval: value < threshold (smaller = more significant).
    For mlog10p: value > threshold (larger = more significant).
    """
    return value > threshold if pval_is_mlog10p else value < threshold

def _passes_threshold_eq(value:float, threshold:float, pval_is_mlog10p:bool)->bool:
    """Like _passes_threshold but inclusive (<=  or >=)."""
    return value >= threshold if pval_is_mlog10p else value <= threshold

def _sort_most_significant_last(items, key, pval_is_mlog10p:bool):
    """Sort so most significant variant is last (for stack popping).
    For raw pval: ascending sort reversed (smallest pval last).
    For mlog10p: ascending sort (largest mlog10p last).
    """
    return sorted(items, key=key, reverse=not pval_is_mlog10p)

def _load_variants_by_chrom(summstat_resource: TabixResource, variants, data_columns:List[str],
                            max_gap:int=1_000_000) -> Dict[Variant,List[str]]:
    """Load data_columns for many variants with one tabix fetch per cluster of nearby
    variants instead of one fetch per variant. Variants on a chromosome are sorted and split
    into clusters wherever the gap between consecutive positions exceeds max_gap, bounding
    each fetch span so a chromosome-wide spread doesn't pull the whole chromosome into memory.
    Returns {Variant: [cols...]} for the requested variants found in the file.
    """
    by_chrom: Dict[str,List[Variant]] = {}
    for v in variants:
        by_chrom.setdefault(v.chrom, []).append(v)
    out: Dict[Variant,List[str]] = {}
    for chrom, vs in by_chrom.items():
        want = set(vs)
        positions = sorted(set(v.pos for v in vs))
        clusters = []
        start = prev = positions[0]
        for p in positions[1:]:
            if p - prev > max_gap:
                clusters.append((start, prev))
                start = p
            prev = p
        clusters.append((start, prev))
        for lo, hi in clusters:
            region = summstat_resource.load_region(chrom, lo, hi+1, data_columns)
            for var, data in region.items():
                if var in want:
                    out[var] = data
    return out

def credible_grouping(credsets:List[CS],summstat_resource: TabixResource,ld_api: LDAccess, options: GroupingOptions) -> Dict[Variant,List[Var]]:
    """Return ld partners for credible sets. Returned as a dictionary of lead var -> ld partner keys & values
    """
    output = {}
    #create list of credible set lead variants
    cs_lead_vars = [a.lead for a in credsets]
    #load pvals from summstat (one fetch per cluster of nearby leads instead of one per lead)
    fetched = _load_variants_by_chrom(summstat_resource, cs_lead_vars, [options.column_names.pval])
    cs_pvals = {}
    for c in cs_lead_vars:
        try:
            cs_pvals[c] = float(fetched[c][0])
        except:
            raise Exception(f"Fatal Error: Credible set lead variant {c} was not found in the summary statistic resource. Can not continue.")
    csleads_by_pval = sorted([(key,value)for key,value in cs_pvals.items()],key=lambda x:x[1],
        reverse=options.pval_is_mlog10p)
    variants_that_can_not_be_included:set[Variant] = set()
    #prefetch LD for all cs leads in parallel (threshold 0), apply dynamic threshold in-memory below
    total_cs = len(csleads_by_pval)
    print(f"credible_grouping: prefetching LD for {total_cs} credible set leads with {options.ld_workers} worker(s)", flush=True)
    ld_cache = ld_api.get_ranges([lv for lv,_ in csleads_by_pval], options.range, workers=options.ld_workers)
    for cs_idx,(lead_var,lead_pval) in enumerate(csleads_by_pval,1):
        print(f"credible_grouping: credible set {cs_idx}/{total_cs}", flush=True)
        # get LD partners
        #dynamic/static r2 thresh decision
        ld_thresh = ld_threshold(options.r2_threshold,options.ld_mode,lead_pval,options.pval_is_mlog10p)
        ld_partners = [a for a in ld_cache[lead_var] if a.r2 > ld_thresh]
        ld_partners = [a for a in ld_partners if a.variant1 == lead_var]
        # get locus in summary statistic
        summstat_data = {}
        for var, vardata in  summstat_resource.load_region(lead_var.chrom,lead_var.pos-options.range,lead_var.pos+options.range,[options.column_names.pval,options.column_names.beta]).items():
            try:
                summstat_data[var] = {"pval":float(vardata[0]),"beta":float(vardata[1])}
            except:
                continue
        #filter by p-value
        summstat_data = {a:b for a,b in summstat_data.items() if _passes_threshold_eq(b["pval"], options.p2_threshold, options.pval_is_mlog10p)}
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
    p1_piles:Dict[str,List[Var]] = {}
    p2_piles:Dict[str,List[Var]] = {}
    for s in summstat_resource.sequences:
        p1_piles[s] = []
        p2_piles[s] = []
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
    # hoist header->index lookups to int locals once instead of a dict lookup per column per row
    i_c, i_p, i_r, i_a, i_pval, i_beta = (hdi[cpra[0]], hdi[cpra[1]], hdi[cpra[2]],
                                          hdi[cpra[3]], hdi[cpra[4]], hdi[cpra[5]])
    for cols in summstat_resource.fetch_all_tuples():
        try:
            pval = float(cols[i_pval])
        except:
            continue
        if _passes_threshold(pval, options.p2_threshold, options.pval_is_mlog10p):
            # beta only needed for rows that pass p2; most genome-wide rows don't
            try:
                beta = float(cols[i_beta])
            except:
                continue
            chrom = cols[i_c]
            var = Var(Variant(
                chrom,
                int(cols[i_p]),
                cols[i_r],
                cols[i_a]),
                pval,
                beta,
                None
            )
            p2_piles[chrom].append(var)
            if _passes_threshold(pval, options.p1_threshold, options.pval_is_mlog10p):
                p1_piles[chrom].append(var)
    for chrom in p1_piles.keys():

        p1_pile = _sort_most_significant_last(p1_piles[chrom],key=lambda x: (x.pval,x.id.chrom,x.id.pos), pval_is_mlog10p=options.pval_is_mlog10p)
        p2_pile = p2_piles[chrom]
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
    return output

def _greedy_ld_group(p1_piles: Dict[str,List[Var]], p2_piles: Dict[str,List[Var]],
                     ld_cache: Dict[Variant,List], options: GroupingOptions) -> List[PeakLocus]:
    """Greedy LD grouping over prefetched LD (cache: lead Variant -> list of LDData).

    Equivalent to the original per-iteration list-rebuild loop, but avoids its
    O(loci * pile size) rescans: p2 is indexed by variant id once, and consumed variants
    are tracked in a set instead of rebuilding the p1/p2 lists on every iteration. Output
    (leads, partner sets, r2 values, partner order) is identical to the old loop.
    """
    output: List[PeakLocus] = []
    total_p1 = sum(len(v) for v in p1_piles.values())
    locus_num = 0
    for seq in p1_piles.keys():
        # most significant variant last, so reversed() walks most -> least significant
        p1_pile = _sort_most_significant_last(p1_piles[seq],key=lambda x: x.pval,
                                              pval_is_mlog10p=options.pval_is_mlog10p)
        p2_pile = p2_piles[seq]
        # index p2 by variant id, keeping original (file) order so partner lists match
        # the old p2-order output exactly
        p2_by_id: Dict[Variant,Var] = {}
        p2_order: Dict[Variant,int] = {}
        for i,a in enumerate(p2_pile):
            p2_by_id[a.id] = a
            p2_order[a.id] = i
        # variants claimed as LD partners of an already-processed (more significant) lead:
        # excluded from being future leads (always) and from being partners again (non-overlap)
        consumed: set[Variant] = set()
        for lead_var in reversed(p1_pile):
            lead_variant = lead_var.id
            if lead_variant in consumed:
                continue
            locus_num += 1
            print(f"ld_grouping: locus {locus_num}/{total_p1} (chr{seq})", flush=True)
            ld_thresh = ld_threshold(options.r2_threshold,options.ld_mode,lead_var.pval,options.pval_is_mlog10p)
            # apply per-lead threshold to the prefetched (threshold-0) cache entry; dict keeps
            # the last r2 per variant2, matching the old index-dict behavior
            r2_by_v2 = {a.variant2:a.r2 for a in ld_cache[lead_variant]
                        if a.r2 > ld_thresh and a.variant1 == lead_variant and a.variant2 != lead_variant}
            # partners must be in p2 and, unless overlap, not already consumed; keep p2 order
            partner_ids = [v2 for v2 in r2_by_v2
                           if v2 in p2_by_id and (options.overlap or v2 not in consumed)]
            partner_ids.sort(key=lambda v: p2_order[v])
            ld_partners = [Var(v2,p2_by_id[v2].pval,p2_by_id[v2].beta,r2_by_v2[v2]) for v2 in partner_ids]
            # every LD partner of this lead is removed from future lead/partner consideration
            consumed.update(r2_by_v2.keys())
            output.append(PeakLocus(lead_var,ld_partners,Grouping.LD))
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
    p1_piles:Dict[str,List[Var]] = {}
    p2_piles:Dict[str,List[Var]] = {}
    for s in summstat_resource.sequences:
        p1_piles[s] = []
        p2_piles[s] = []
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
    # hoist header->index lookups to int locals once instead of a dict lookup per column per row
    i_c, i_p, i_r, i_a, i_pval, i_beta = (hdi[cpra[0]], hdi[cpra[1]], hdi[cpra[2]],
                                          hdi[cpra[3]], hdi[cpra[4]], hdi[cpra[5]])
    for cols in summstat_resource.fetch_all_tuples():
        try:
            pval = float(cols[i_pval])
        except:
            continue
        if _passes_threshold(pval, options.p2_threshold, options.pval_is_mlog10p):
            # beta only needed for rows that pass p2; most genome-wide rows don't
            try:
                beta = float(cols[i_beta])
            except:
                continue
            chrom = cols[i_c]
            var = Var(Variant(
                chrom,
                int(cols[i_p]),
                cols[i_r],
                cols[i_a]),
                pval,
                beta,
                1.0
            )
            p2_piles[chrom].append(var)
            if _passes_threshold(pval, options.p1_threshold, options.pval_is_mlog10p):
                p1_piles[chrom].append(var)
    # prefetch LD for every candidate lead up front (in parallel); the greedy loop below then
    # reads from cache rather than blocking on one LD fetch per peak. LD for a lead is
    # independent of grouping state, so this is safe. Fetched at threshold 0 so the per-lead
    # dynamic threshold can be applied in memory.
    total_p1 = sum(len(v) for v in p1_piles.values())
    all_leads = [v.id for seq in p1_piles for v in p1_piles[seq]]
    print(f"ld_grouping: prefetching LD for {total_p1} candidate lead variants with {options.ld_workers} worker(s)", flush=True)
    ld_cache = ld_api.get_ranges(all_leads, options.range, workers=options.ld_workers)
    return _greedy_ld_group(p1_piles, p2_piles, ld_cache, options)

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
    # hoist header->index lookups to int locals once instead of a dict lookup per column per row
    i_c, i_p, i_r, i_a, i_pval, i_beta = (hdi[cpra[0]], hdi[cpra[1]], hdi[cpra[2]],
                                          hdi[cpra[3]], hdi[cpra[4]], hdi[cpra[5]])
    for cols in summstat_resource.fetch_all_tuples():
        try:
            pval = float(cols[i_pval])
        except:
            continue
        if _passes_threshold(pval, options.p1_threshold, options.pval_is_mlog10p):
            # beta only needed for rows that pass the threshold
            try:
                beta = float(cols[i_beta])
            except:
                continue
            var = Var(Variant(
                cols[i_c],
                int(cols[i_p]),
                cols[i_r],
                cols[i_a]),
                pval,
                beta,
                None
            )
            p1_pile.append(var)
    return p1_pile

@timefunc
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
        #one fetch per cluster of nearby cs variants instead of one per variant
        cs_pvalbetadict = {}
        fetched_cs = _load_variants_by_chrom(summstat_resource, cs_variants,
            [gr_opts.column_names.pval,gr_opts.column_names.beta])
        for var, data in fetched_cs.items():
            try:
                cs_pvalbetadict[var] = [float(d) for d in data]
            except:
                continue
        cs_groups = {}
        cs_infos = {}
        #prefetch internal-CS LD per lead (each cs uses its own range), in parallel
        lead_ranges = {}
        for c in cs:
            rng = max([abs(c.lead.pos - a.variant.pos) for a in c.variants])+5000
            lead_ranges[c.lead] = max(lead_ranges.get(c.lead,0), rng)
        cs_lead_list = list(lead_ranges.keys())
        print(f"cs_grouping: prefetching LD for {len(cs_lead_list)} credible set leads with {gr_opts.ld_workers} worker(s)", flush=True)
        cs_ld_cache = ld_access.get_ranges(cs_lead_list, 0,
            bp_ranges=[lead_ranges[l] for l in cs_lead_list], workers=gr_opts.ld_workers)
        #get pval, beta, r2 for credible set variants
        for c in cs:
            #gather ld data
            c_variants = set([a.variant for a in c.variants])
            #susie might not have LD so it was prefetched above
            ld_data = cs_ld_cache[c.lead]
            ld_data = [a for a in ld_data if ( (a.variant1 == c.lead) and (a.variant2 in c_variants))  ]
            ld_dict = {a.variant2:a.r2 for a in ld_data}

            vs = [Var(
                    v.variant,
                    cs_pvalbetadict.get(v.variant,(float("nan"),float("nan")))[0],
                    cs_pvalbetadict.get(v.variant,(float("nan"),float("nan")))[1],
                    ld_dict.get(v.variant,float("nan"))
                )
                for v in c.variants if v.variant != c.lead]
            lead_var = Var(c.lead,cs_pvalbetadict.get(c.lead,(float("nan"),float("nan")))[0],
                    cs_pvalbetadict.get(c.lead,(float("nan"),float("nan")))[1],
                    1.0)
            #add lead var to cs vars
            vs.append(lead_var)
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
