#include "PersonBuilder.h"
#include "PersonAddressBuilder.h"
#include "PersonJobBuilder.h"

PersonBuilderBase::PersonBuilderBase(Person& person) : person(person) {}

// In the below constructor, we are initializing the variable in base class of person builder with its own variable. Base class takes a reference.
PersonBuilder::PersonBuilder() :PersonBuilderBase(p) {};

PersonAddressBuilder PersonBuilderBase::lives() const
{
	return PersonAddressBuilder{ person};
}

PersonJobBuilder PersonBuilderBase::works() const
{
	return PersonJobBuilder{ person };
}


