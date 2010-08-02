template <typename TSimType>  class DefaultCodFieldTraits
{
   typedef PhysicalScalarBase<TSimType> OCFieldBase;
   typedef PhysicalScalarBase<TSimType> ICFieldBase;

   typedef PhysicalScalar<TSimType, OCFieldBase> OCFieldType;
   typedef PhysicalScalar<TSimType, ICFieldBase> ICFieldType;
}
