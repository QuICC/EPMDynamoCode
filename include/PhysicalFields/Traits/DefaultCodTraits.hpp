template <typename TSimType>  class DefaultCodTraits
{
   typedef DefaultCodDomainTraits<TSimType, FieldTraits> DomainTraits;

   typedef OuterCoreSource<TSimType, DomainTraits> VariableType;
}
