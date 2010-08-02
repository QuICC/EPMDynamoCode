template <typename TSimType>  class DefaultCodDomainTraits
{
   typedef DefaultCodFieldTraits<TSimType>   FieldTraits;

   typedef FullSphereField<TSimType, FieldTraits> DomainType;
}
