#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <string>
using namespace std;

#include "Person.h"
#include "PersonBuilder.h"
#include "PersonAddressBuilder.h"
#include "PersonJobBuilder.h"

int main()
{
    Person p = Person::create()         // Creates a PersonBuilder object, which inturn creates a person object as one of its attributes.
        .lives().at("123 London Road")  // Creates a PersonAddressBuilder object by taking in person object as reference and then filling that object address details.
            .with_postcode("SW1 1GB")
            .in("London")
        .works().at("PragmaSoft")       // Creates a PersonJobBuilder object by taking in person object as reference and then filling that object income details.
            .as_a("Consultant")
            .earning(10e6);

    cout << p << endl;
    getchar();
    return 0;
}