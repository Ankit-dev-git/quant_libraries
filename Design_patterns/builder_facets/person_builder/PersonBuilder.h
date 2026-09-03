#pragma once
#include "Person.h"

class PersonBuilderBase {
protected:
	// This class will have Person object reference as one of its variable. 
	// In the constructor, we will pass the object created from Person class to ensure no duplication of objects created in builder classes,
	// as builder classes inherit from this class and not PersonBuilder class
	Person& person;
public:
	PersonBuilderBase(Person& P);

	operator Person() const {
		return std::move(person);
	}

	PersonAddressBuilder lives() const;
	PersonJobBuilder works() const;
};

class PersonBuilder: public PersonBuilderBase {
	// This class constructor will create the person object. Other classes will just take this object by reference
	Person p;
public:
	PersonBuilder();

};
