/*
 * ShallowRectangle.h  27.09.2006 (mueller)
 *
 * Copyright (C) 2006 by Universitaet Stuttgart (VIS). Alle Rechte vorbehalten.
 */

#pragma once
#if defined(_WIN32) && defined(_MANAGED)
#pragma managed(push, off)
#endif /* defined(_WIN32) && defined(_MANAGED) */


#include "vislib/math/AbstractRectangle.h"


namespace vislib::math {


/**
 * This class represents a shallow rectangle, that uses memory provided
 * by the caller.
 */
template<class T>
class ShallowRectangle : public AbstractRectangle<T, T*> {

public:
    /**
     * Construct a rectangle from an array holding its bounds. The array
     * 'bounds' holds in this order to following borders of the rectangle:
     * left, bottom, right, top.
     *
     * @param bounds The left, bottom, right and top border of the
     *               rectangle in a consecutive order.
     */
    explicit inline ShallowRectangle(T* bounds) {
        this->bounds = bounds;
    }

    /**
     * Copy ctor. This ctor creates an alias!
     *
     * @param rhs The object to clone.
     */
    inline ShallowRectangle(const ShallowRectangle& rhs) {
        this->bounds = rhs.bounds;
    }

    /** Dtor. */
    ~ShallowRectangle();

    /**
     * Assigment operator. This operator never creates an alias, even for
     * shallow rectangles!
     *
     * @param rhs The right hand side operand.
     *
     * @return *this.
     */
    inline ShallowRectangle& operator=(const ShallowRectangle& rhs) {
        Super::operator=(rhs);
        return *this;
    }

    /**
     * Assigment operator. This operator never creates an alias, even for
     * shallow rectangles!
     *
     * This assignment allows for arbitrary rectangle to rectangle
     * conversions.
     *
     * @param rhs The right hand side operand.
     *
     * @return *this.
     */
    template<class Tp, class Sp>
    inline ShallowRectangle& operator=(const AbstractRectangle<Tp, Sp>& rhs) {
        Super::operator=(rhs);
        return *this;
    }

public:
    /** Typedef for the super class. */
    typedef AbstractRectangle<T, T*> Super;
};


/*
 * vislib::math::ShallowRectangle<T>::~ShallowRectangle
 */
template<class T>
ShallowRectangle<T>::~ShallowRectangle() {}


} // namespace vislib::math

#if defined(_WIN32) && defined(_MANAGED)
#pragma managed(pop)
#endif /* defined(_WIN32) && defined(_MANAGED) */
