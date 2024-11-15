from typing import List, Dict, NamedTuple,Any,Optional, Union
from data_access.db import Variant
import abc

#a single annotation for a single variant
Col = str
Value = Union[str,int,float,Variant,None]
VariantAnnotation = Dict[Col,Value]
#a single variant might have multiple instances of annotations, e.g. if there are multiple gwas catalog its for a single variant
Annotation = List[VariantAnnotation]

class AnnotationSource:
    @abc.abstractmethod
    def annotate_variants(self,variants:set[Variant])->Dict[Variant,Annotation]:
        pass

    @staticmethod
    @abc.abstractmethod
    def get_name() -> str:
        """Unique identifier for the annotation
        """
        pass

    @staticmethod
    @abc.abstractmethod
    def get_output_columns() -> List[str]:
        pass