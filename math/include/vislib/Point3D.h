/*
 * Point2D.h
 *
 * Copyright (C) 2006 by Universitaet Stuttgart (VIS). Alle Rechte vorbehalten.
 * Copyright (C) 2005 by Christoph Mueller. Alle Rechte vorbehalten.
 */

#ifndef VISLIB_POINT3D_H_INCLUDED
#define VISLIB_POINT3D_H_INCLUDED
#if _MSC_VER > 1000
#pragma once
#endif /* _MSC_VER > 1000 */


#include "vislib/AbstractPoint3D.h"
#include "vislib/Point.h"


namespace vislib {
namespace math {

    /**
     * A point in three-dimensional space.
     *
     * @author Christoph Mueller
     */
    template<class T >
    class Point3D : public AbstractPoint3D<T, T[3]> {

    public:

        /** A typedef for the super class having the storage related info. */
        typedef AbstractPoint3D<T, T[3]> Super;

        /** 
         * Create a point in the origin of the coordinate system. 
         */
        inline Point3D(void) {
            this->coordinates[0] = static_cast<T>(0);
            this->coordinates[1] = static_cast<T>(0);
            this->coordinates[2] = static_cast<T>(0);
        }

        /**
         * Create a new point.
         *
         * @param x The x-coordinate.
         * @param y The y-coordinate.
         * @param z The z-coordinate.
         */
        inline Point3D(const T& x, const T& y, const T& z) {
            this->coordinates[0] = x;
            this->coordinates[1] = y;
            this->coordinates[2] = z;
        }

        /**
         * Create a new point with the specified coordinates. 'coordinates' must be
         * an array with at least three elements. The caller remains owner of
         * the memory, a deep copy is created.
         *
         * @param coordinates The coordinates of the point.
         */
        inline explicit Point3D(const T *coordinates) {
            ASSERT(coordinates != NULL);
            ::memcpy(this->coordinates, coordinates, 3 * sizeof(T));
        }

        /**
         * Create a new point from its position vector.
         *
         * @param posVec The position vector of the point.
         */
        template<class Tp, class Sp>
        inline explicit Point3D(const AbstractVector<Tp, 3, Sp>& posVec) {
            this->coordinates[0] = posVec[0];
            this->coordinates[1] = posVec[1];
        }

        /**
         * Clone 'rhs'.
         *
         * @param rhs The object to be cloned.
         */
        inline Point3D(const Point3D& rhs) {
            ::memcpy(this->coordinates, rhs.coordinates, 3 * sizeof(T));
        }

        /**
         * Create a copy of 'rhs'. This ctor allows for arbitrary point to
         * point conversions. See documentation of Point for additional
         * remarks.
         *
         * @param rhs The object to be cloned.
         */        
        template<class Tp, unsigned int Dp, class Sp>
        inline Point3D(const AbstractPoint<Tp, Dp, Sp>& rhs) {
            this->coordinates[0] = (Dp < 1) ? static_cast<T>(0) 
                                            : static_cast<T>(rhs[0]);
            this->coordinates[1] = (Dp < 2) ? static_cast<T>(0) 
                                            : static_cast<T>(rhs[1]);
            this->coordinates[2] = (Dp < 3) ? static_cast<T>(0) 
                                            : static_cast<T>(rhs[2]);
        }

        /** Dtor. */
        ~Point3D(void);


        /**
         * Assignment operator.
         *
         * @param rhs The right hand side operand.
         *
         * @return *this.
         */
        inline Point3D& operator =(const Point3D& rhs) {
            Super::operator =(rhs);
            return *this;
        }

        /**
         * Assigment for arbitrary points. A valid static_cast between T and Tp
         * is a precondition for instantiating this template.
         *
         * This operation does <b>not</b> create aliases. 
         *
         * If the two operands have different dimensions, the behaviour is as 
         * follows: If the left hand side operand has lower dimension, the 
         * highest (Dp - D) dimensions are discarded. If the left hand side
         * operand has higher dimension, the missing dimensions are filled with 
         * zero components.
         *
         * Subclasses must ensure that sufficient memory for the 'coordinates'
         * member has been allocated before calling this operator.
         *
         * @param rhs The right hand side operand.
         *
         * @return *this
         */
        template<class Tp, unsigned int Dp, class Sp>
        inline Point3D& operator =(const AbstractPoint<Tp, Dp, Sp>& rhs) {
            Super::operator =(rhs);
            return *this;
        }
    };


    /*
     * vislib::math::Point3D<T>::~Point3D
     */
    template<class T>
    Point3D<T>::~Point3D(void) {
    }

} /* end namespace math */
} /* end namespace vislib */

#endif /* VISLIB_POINT3D_H_INCLUDED */
