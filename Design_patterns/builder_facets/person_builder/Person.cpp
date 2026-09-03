#include "Person.h"
#include "PersonBuilder.h" // Add this include to resolve incomplete type

PersonBuilder Person::create()
{
    return PersonBuilder();
}
