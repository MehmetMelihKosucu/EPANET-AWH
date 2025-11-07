/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

//! \file element.h
//! \brief Describes the Element class.

#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <string>

//! \class Element
//! \brief Abstract parent class for all pipe network components.

class Element
{
  public:

    enum ElementType {NODE, LINK, PATTERN, CURVE, CONTROL};

    Element(std::string name_);
    virtual ~Element() = 0;

    std::string name;       //!< element's ID name
    int         index;      //!< index in array of elements

    // Add a public copy constructor
    Element(const Element& e) = default;

  private:

    // Elements can't be copied or tested for equality
    
    Element& operator=(const Element& e);
};

#endif
