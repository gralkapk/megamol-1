/*
 * AbstractSimpleMessageHeader.h
 *
 * Copyright (C) 2006 - 2010 by Visualisierungsinstitut Universitaet Stuttgart.
 * Alle Rechte vorbehalten.
 */

#pragma once
#if defined(_WIN32) && defined(_MANAGED)
#pragma managed(push, off)
#endif /* defined(_WIN32) && defined(_MANAGED) */


#include "vislib/net/SimpleMessageHeaderData.h"


namespace vislib::net {


/**
 * This class provides the interface for a message header in the simple
 * network protocol implementation of VISlib. There are two subclasses of
 * this abstract class: SimpleMessageHeader is an implementation which
 * provides storage for the header data. ShallowSimpleMessageHeader does
 * not provide any storage, but takes a pointer to memory containing the
 * header data.
 */
class AbstractSimpleMessageHeader {

public:
    /** Dtor. */
    virtual ~AbstractSimpleMessageHeader();

    /**
     * Answer the body size stored in the message header.
     *
     * @return The body size.
     */
    inline SimpleMessageSize GetBodySize() const {
        return this->PeekData()->BodySize;
    }

    /**
     * Answer the size of the header packet. This is the size of the data
     * returned by PeekData().
     *
     * @return The size of the header data in bytes.
     */
    inline SimpleMessageSize GetHeaderSize() const {
        return sizeof(SimpleMessageHeaderData);
    }

    /**
     * Answer the message ID.
     *
     * @return The message ID.
     */
    inline SimpleMessageID GetMessageID() const {
        return this->PeekData()->MessageID;
    }

    /**
     * Answer whether the body size is not zero.
     *
     * @return true if the body size is larger than zero, false otherwise.
     */
    inline bool HasBody() const {
        return (this->PeekData()->BodySize > 0);
    }

    /**
     * Provides direct access to the underlying SimpleMessageHeaderData.
     *
     * @return A pointer to the message header data.
     */
    virtual SimpleMessageHeaderData* PeekData() = 0;

    /**
     * Provides direct access to the underlying SimpleMessageHeaderData.
     *
     * @return A pointer to the message header data.
     */
    virtual const SimpleMessageHeaderData* PeekData() const = 0;

    /**
     * Set the body size.
     *
     * @param bodySize The body size.
     */
    inline void SetBodySize(const SimpleMessageSize bodySize) {
        this->PeekData()->BodySize = bodySize;
    }

    /**
     * Set a new message ID.
     *
     * @param messageID  The new message ID.
     * @param isSystemID Disables the system ID check. Must be false.
     *
     * @throw IllegalParamException If the message ID is a system ID.
     */
    void SetMessageID(const SimpleMessageID messageID, bool isSystemID = false);

    /**
     * Assignment operator.
     *
     * @param The right hand side operand.
     *
     * @return *this
     */
    AbstractSimpleMessageHeader& operator=(const AbstractSimpleMessageHeader& rhs);

    /**
     * Assignment operator.
     *
     * @param The right hand side operand.
     *
     * @return *this
     */
    AbstractSimpleMessageHeader& operator=(const SimpleMessageHeaderData& rhs);

    /**
     * Assignment operator.
     *
     * @param The right hand side operand. This must not be NULL.
     *
     * @return *this
     */
    AbstractSimpleMessageHeader& operator=(const SimpleMessageHeaderData* rhs);

    /**
     * Test for equality.
     *
     * @param The right hand side operand.
     *
     * @return true in case this object and 'rhs' are equal,
     *         false otherwise.
     */
    bool operator==(const AbstractSimpleMessageHeader& rhs) const;

    /**
     * Test for inequality.
     *
     * @param The right hand side operand.
     *
     * @return true in case this object and 'rhs' are not equal,
     *         false otherwise.
     */
    inline bool operator!=(const AbstractSimpleMessageHeader& rhs) const {
        return !(*this == rhs);
    }

protected:
    /** Ctor. */
    AbstractSimpleMessageHeader();
};

} // namespace vislib::net

#if defined(_WIN32) && defined(_MANAGED)
#pragma managed(pop)
#endif /* defined(_WIN32) && defined(_MANAGED) */
