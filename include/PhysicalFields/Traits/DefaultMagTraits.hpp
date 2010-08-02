template <typename TSimType>  class DefaultMagTraits
{
   typedef  DefaultMagFieldTraits<TSimType>   FieldTraits;

   typedef FullSphereField<TSimType, FieldTraits> VariableType;
}
