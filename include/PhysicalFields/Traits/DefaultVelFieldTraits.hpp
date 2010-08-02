template <typename TSimType>  class DefaultVelFieldTraits
{
   typedef PhysicalTorPolBase<TSimType> OCFieldBase;
   typedef PhysicalTorPolBase<TSimType> ICFieldBase;

   typedef PhysicalTorPolField<TSimType, OCFieldBase> OCFieldType;
   typedef PhysicalTorPolField<TSimType, ICFieldBase> ICFieldType;
}
