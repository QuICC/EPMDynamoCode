template <typename TSimType>  class DefaultVelTraits
{
   typedef  DefaultVelFieldTraits<TSimType>   FieldTraits;

   typedef FullSphereField<TSimType, FieldTraits> VariableType;
}
