/*
 * MMSPDDataSource.cpp
 *
 * Copyright (C) 2011 by VISUS (Universitaet Stuttgart)
 * Alle Rechte vorbehalten.
 */

#include "stdafx.h"
#include "MMSPDDataSource.h"
#include "param/FilePathParam.h"
#include "MultiParticleDataCall.h"
#include "CoreInstance.h"
//#include "vislib/Log.h"
#include "vislib/File.h"
#include "vislib/String.h"
#include "vislib/StringTokeniser.h"
//#include "vislib/SystemInformation.h"
#include "vislib/mathfunctions.h"
#include "vislib/UTF8Encoder.h"
#include "vislib/VersionNumber.h"

using namespace megamol::core;


/* defines for the frame cache size */
// minimum number of frames in the cache (2 for interpolation; 1 for loading)
#define CACHE_SIZE_MIN 3
// maximum number of frames in the cache (just a nice number)
#define CACHE_SIZE_MAX 1000
// factor multiplied to the frame size for estimating the overhead to the pure data.
#define CACHE_FRAME_FACTOR 1.2f

/*****************************************************************************/

/*
 * moldyn::MMSPDDataSource::Frame::Frame
 */
moldyn::MMSPDDataSource::Frame::Frame(view::AnimDataModule& owner)
        : moldyn::MMSPDFrameData(), view::AnimDataModule::Frame(owner) {
    // intentionally empty
}


/*
 * moldyn::MMSPDDataSource::Frame::~Frame
 */
moldyn::MMSPDDataSource::Frame::~Frame() {
    // intentionally empty
}


///*
// * moldyn::MMPLDDataSource::Frame::LoadFrame
// */
//bool moldyn::MMPLDDataSource::Frame::LoadFrame(vislib::sys::File *file, unsigned int idx, UINT64 size) {
//    this->frame = idx;
//    this->dat.EnforceSize(static_cast<SIZE_T>(size));
//    return (file->Read(this->dat, size) == size);
//}
//
//
///*
// * moldyn::MMPLDDataSource::Frame::SetData
// */
//void moldyn::MMPLDDataSource::Frame::SetData(MultiParticleDataCall& call) {
//    if (this->dat.IsEmpty()) {
//        call.SetParticleListCount(0);
//        return;
//    }
//
//    SIZE_T p = sizeof(UINT32);
//    UINT32 plc = *this->dat.As<UINT32>();
//    call.SetParticleListCount(plc);
//    for (UINT32 i = 0; i < plc; i++) {
//        MultiParticleDataCall::Particles &pts = call.AccessParticles(i);
//
//        UINT8 vrtType = *this->dat.AsAt<UINT8>(p); p += 1;
//        UINT8 colType = *this->dat.AsAt<UINT8>(p); p += 1;
//        MultiParticleDataCall::Particles::VertexDataType vrtDatType;
//        MultiParticleDataCall::Particles::ColourDataType colDatType;
//        SIZE_T vrtSize = 0;
//        SIZE_T colSize = 0;
//
//        switch (vrtType) {
//            case 0: vrtSize = 0; vrtDatType = MultiParticleDataCall::Particles::VERTDATA_NONE; break;
//            case 1: vrtSize = 12; vrtDatType = MultiParticleDataCall::Particles::VERTDATA_FLOAT_XYZ; break;
//            case 2: vrtSize = 16; vrtDatType = MultiParticleDataCall::Particles::VERTDATA_FLOAT_XYZR; break;
//            case 3: vrtSize = 6; vrtDatType = MultiParticleDataCall::Particles::VERTDATA_SHORT_XYZ; break;
//            default: vrtSize = 0; vrtDatType = MultiParticleDataCall::Particles::VERTDATA_NONE; break;
//        }
//        if (vrtType != 0) {
//            switch (colType) {
//                case 0: colSize = 0; colDatType = MultiParticleDataCall::Particles::COLDATA_NONE; break;
//                case 1: colSize = 3; colDatType = MultiParticleDataCall::Particles::COLDATA_UINT8_RGB; break;
//                case 2: colSize = 4; colDatType = MultiParticleDataCall::Particles::COLDATA_UINT8_RGBA; break;
//                case 3: colSize = 4; colDatType = MultiParticleDataCall::Particles::COLDATA_FLOAT_I; break;
//                case 4: colSize = 12; colDatType = MultiParticleDataCall::Particles::COLDATA_FLOAT_RGB; break;
//                case 5: colSize = 16; colDatType = MultiParticleDataCall::Particles::COLDATA_FLOAT_RGBA; break;
//                default: colSize = 0; colDatType = MultiParticleDataCall::Particles::COLDATA_NONE; break;
//            }
//        } else {
//            colDatType = MultiParticleDataCall::Particles::COLDATA_NONE;
//            colSize = 0;
//        }
//        unsigned int stride = static_cast<unsigned int>(vrtSize + colSize);
//
//        if ((vrtType == 1) || (vrtType == 3)) {
//            pts.SetGlobalRadius(*this->dat.AsAt<float>(p)); p += 4;
//        } else {
//            pts.SetGlobalRadius(0.05f);
//        }
//
//        if (colType == 0) {
//            pts.SetGlobalColour(*this->dat.AsAt<UINT8>(p),
//                *this->dat.AsAt<UINT8>(p + 1),
//                *this->dat.AsAt<UINT8>(p + 2));
//            p += 4;
//        } else {
//            pts.SetGlobalColour(192, 192, 192);
//            if (colType == 3) {
//                pts.SetColourMapIndexValues(
//                    *this->dat.AsAt<float>(p),
//                    *this->dat.AsAt<float>(p + 4));
//                p += 8;
//            } else {
//                pts.SetColourMapIndexValues(0.0f, 1.0f);
//            }
//        }
//
//        pts.SetCount(*this->dat.AsAt<UINT64>(p)); p += 8;
//
//        pts.SetVertexData(vrtDatType, this->dat.At(p), stride);
//        pts.SetColourData(colDatType, this->dat.At(p + vrtSize), stride);
//
//    }
//
//}

/*****************************************************************************/

/*
 * moldyn::MMSPDDataSource::FileFormatAutoDetect
 */
float moldyn::MMSPDDataSource::FileFormatAutoDetect(const unsigned char* data, SIZE_T dataSize) {
    return (((dataSize >= 6)
        && ((::memcmp(data, "MMSPDb", 6) == 0)
            || (::memcmp(data, "MMSPDa", 6) == 0)
            || (::memcmp(data, "MMSPDu", 6) == 0)))
        || ((dataSize >= 9)
        && (::memcmp(data, "\xEF\xBB\xBFMMSPDu", 9) == 0))) ? 1.0f : 0.0f;
}


/*
 * moldyn::MMSPDDataSource::MMSPDDataSource
 */
moldyn::MMSPDDataSource::MMSPDDataSource(void) : view::AnimDataModule(),
        filename("filename", "The path to the MMSPD file to load."),
        getData("getdata", "Slot to request data from this data source."),
        dataHeader(), file(NULL), frameIdx(NULL), framePartCnts(NULL),
        clipbox(-1.0f, -1.0f, -1.0f, 1.0f, 1.0f, 1.0f) {

    this->filename.SetParameter(new param::FilePathParam(""));
    this->filename.SetUpdateCallback(&MMSPDDataSource::filenameChanged);
    this->MakeSlotAvailable(&this->filename);

    this->getData.SetCallback("MultiParticleDataCall", "GetData", &MMSPDDataSource::getDataCallback);
    this->getData.SetCallback("MultiParticleDataCall", "GetExtent", &MMSPDDataSource::getExtentCallback);
    this->MakeSlotAvailable(&this->getData);

    this->setFrameCount(1);
    this->initFrameCache(1);
}


/*
 * moldyn::MMSPDDataSource::~MMSPDDataSource
 */
moldyn::MMSPDDataSource::~MMSPDDataSource(void) {
    this->Release();
}


/*
 * moldyn::MMSPDDataSource::constructFrame
 */
view::AnimDataModule::Frame* moldyn::MMSPDDataSource::constructFrame(void) const {
    Frame *f = new Frame(*const_cast<moldyn::MMSPDDataSource*>(this));
    return f;
}


/*
 * moldyn::MMSPDDataSource::create
 */
bool moldyn::MMSPDDataSource::create(void) {
    return true;
}


/*
 * moldyn::MMSPDDataSource::loadFrame
 */
void moldyn::MMSPDDataSource::loadFrame(view::AnimDataModule::Frame *frame,
        unsigned int idx) {
    using vislib::sys::Log;
    Frame *f = dynamic_cast<Frame*>(frame);
    if (f == NULL) return;
    if (this->file == NULL) {
        //f->Clear();
        return;
    }
    //printf("Requesting frame %u of %u frames\n", idx, this->FrameCount());
    ASSERT(idx < this->FrameCount());
    this->file->Seek(this->frameIdx[idx]);
    //if (!f->LoadFrame(this->file, idx, this->frameIdx[idx + 1] - this->frameIdx[idx])) {
    //    // failed
    //    Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, "Unable to read frame %d from MMPLD file\n", idx);
    //}
}


/*
 * moldyn::MMSPDDataSource::release
 */
void moldyn::MMSPDDataSource::release(void) {
    this->resetFrameCache();
    if (this->file != NULL) {
        vislib::sys::File *f = this->file;
        this->file = NULL;
        f->Close();
        delete f;
    }
    ARY_SAFE_DELETE(this->frameIdx);
    ARY_SAFE_DELETE(this->framePartCnts);
}



/*
 * moldyn::MMSPDDataSource::filenameChanged
 */
bool moldyn::MMSPDDataSource::filenameChanged(param::ParamSlot& slot) {
    using vislib::sys::Log;
    using vislib::sys::File;
    this->resetFrameCache();
    this->dataHeader.SetParticleCount(0);
    this->dataHeader.SetTimeCount(1);
    this->dataHeader.BoundingBox().Set(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0);
    this->dataHeader.Types().Clear();
    this->clipbox = this->dataHeader.GetBoundingBox();
    ARY_SAFE_DELETE(this->frameIdx);
    ARY_SAFE_DELETE(this->framePartCnts);

    if (this->file == NULL) {
        this->file = new vislib::sys::File();
    } else {
        this->file->Close();
    }
    ASSERT(this->filename.Param<param::FilePathParam>() != NULL);

    if (!this->file->Open(this->filename.Param<param::FilePathParam>()->Value(), File::READ_ONLY, File::SHARE_READ, File::OPEN_ONLY)) {
        this->GetCoreInstance()->Log().WriteMsg(Log::LEVEL_ERROR, "Unable to open MMSPD-File \"%s\".", vislib::StringA(
            this->filename.Param<param::FilePathParam>()->Value()).PeekBuffer());

        SAFE_DELETE(this->file);
        this->setFrameCount(1);
        this->initFrameCache(1);

        return true;
    }

#define _ERROR_OUT(MSG) { Log::DefaultLog.WriteMsg(Log::LEVEL_ERROR, MSG); \
        SAFE_DELETE(this->file); \
        this->setFrameCount(1); \
        this->initFrameCache(1); \
        this->dataHeader.SetParticleCount(0); \
        this->dataHeader.SetTimeCount(1); \
        this->dataHeader.BoundingBox().Set(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0); \
        this->dataHeader.Types().Clear(); \
        this->clipbox = this->dataHeader.GetBoundingBox(); \
        ARY_SAFE_DELETE(this->frameIdx); \
        ARY_SAFE_DELETE(this->framePartCnts); \
        return true; }
#define _ASSERT_READFILE(BUFFER, BUFFERSIZE) { if (this->file->Read((BUFFER), (BUFFERSIZE)) != (BUFFERSIZE)) { \
        _ERROR_OUT("Unable to read MMSPD file: seems truncated"); \
    } }
#define _ASSERT_READSTRINGBINARY(STRING) { (STRING).Clear(); while(true) { char c; _ASSERT_READFILE(&c, 1) if (c == 0) break; (STRING).Append(c); } }
#define _ASSERT_READLINE(STRING) { (STRING).Clear(); while(true) { char c; _ASSERT_READFILE(&c, 1) (STRING).Append(c); if (c == 0x0A) break; } }

    // reading format marker
    BYTE headerID[9];
    _ASSERT_READFILE(headerID, 9);
    bool jmpBk, text, unicode, bigEndian;
    if ((text = (::memcmp(headerID, "MMSPDb", 6) != 0))
            && (unicode = (::memcmp(headerID, "MMSPDa", 6) != 0))
            && (::memcmp(headerID, "MMSPDu", 6) != 0)
            && (jmpBk = (::memcmp(headerID, "\xEF\xBB\xBFMMSPDu", 9) != 0))) {
        _ERROR_OUT("MMSPD format marker not found");
    }
    if (jmpBk) {
        this->file->Seek(-3, vislib::sys::File::CURRENT);
    }

    // read format version information
    vislib::VersionNumber version;
    if (text) {
        // the is ultimatively slow, but, it is okey here
        char c;
        vislib::StringA verStr;
        do {
            _ASSERT_READFILE(&c, 1);
            verStr.Append(c);
        } while (c != 0x0A);
        verStr.TrimSpaces();
        version.Parse(verStr);

        if (version != vislib::VersionNumber(1, 0)) {
            vislib::StringA msg;
            msg.Format("Version %s found. Supporting only version 1.0", version.ToStringA().PeekBuffer());
            _ERROR_OUT(msg);
        }

    } else {
        unsigned char buf[14];
        _ASSERT_READFILE(&buf, 14);
        if ((buf[0] != 0x00) || (buf[1] != 0xFF)) {
            _ERROR_OUT("MMSPD file format marker sequence broken @1");
        }
        unsigned int endianessTest;
        ::memcpy(&endianessTest, &buf[2], 4);
        if (endianessTest == 0x78563412) { // which is 2018915346u as specified
            bigEndian = false;
        } else if (endianessTest == 0x12345678) {
            bigEndian = true;
        } else {
            _ERROR_OUT("MMSPD file format marker sequence broken @2");
        }
        unsigned short majorVer, minorVer;
        if (bigEndian) {
            vislib::Swap(buf[6], buf[7]);
            vislib::Swap(buf[8], buf[9]);
        }
        ::memcpy(&majorVer, &buf[6], 2);
        ::memcpy(&minorVer, &buf[8], 2);
        if ((majorVer != 1) && (minorVer != 0)) {
            vislib::StringA msg;
            msg.Format("Version %d.%d found. Supporting only version 1.0", static_cast<int>(majorVer), static_cast<int>(minorVer));
            _ERROR_OUT(msg);
        }
        version.Set(majorVer, minorVer);
        for (int i = 10; i < 14; i++) {
            if (buf[i] < 128) {
                vislib::sys::Log::DefaultLog.WriteWarn("MMSPD file format marker binary guard byte %d illegal", (i - 9));
            }
        }

    }
    // file format marker successfully read
    // file pointer is a the start of the next line/data block

    // reading header line and particle types definitions
    if (text) {
        // reading header line
        vislib::StringA line;
        _ASSERT_READLINE(line);
        if (unicode) {
            vislib::StringW uniLine;
            if (!vislib::UTF8Encoder::Decode(uniLine, line)) {
                _ERROR_OUT("Failed to decode UTF8 header line");
            }
            line = uniLine;
        }
        line.TrimSpaces();
        vislib::Array<vislib::StringA> tokens(vislib::StringTokeniserA::Split(line, ' ', true));
        if (tokens.Count() < 10) {
            _ERROR_OUT("Header line incomplete");
        } else if (tokens.Count() > 10) {
            vislib::sys::Log::DefaultLog.WriteWarn("Trailing information on header line will be ignored");
        }

        const char *fieldName = "unknown";
        try {
            fieldName = "hasIDs";
            bool hasIDs = vislib::CharTraitsA::ParseBool(tokens[0]);
            fieldName = "minX";
            double minX = vislib::CharTraitsA::ParseDouble(tokens[1]);
            fieldName = "minY";
            double minY = vislib::CharTraitsA::ParseDouble(tokens[2]);
            fieldName = "minZ";
            double minZ = vislib::CharTraitsA::ParseDouble(tokens[3]);
            fieldName = "maxX";
            double maxX = vislib::CharTraitsA::ParseDouble(tokens[4]);
            fieldName = "maxY";
            double maxY = vislib::CharTraitsA::ParseDouble(tokens[5]);
            fieldName = "maxZ";
            double maxZ = vislib::CharTraitsA::ParseDouble(tokens[6]);
            fieldName = "timeCount";
            UINT32 timeCount = static_cast<UINT32>(vislib::CharTraitsA::ParseUInt64(tokens[7]));
            fieldName = "typeCount";
            UINT32 typeCount = static_cast<UINT32>(vislib::CharTraitsA::ParseUInt64(tokens[8]));
            fieldName = "partCount";
            UINT64 partCount = vislib::CharTraitsA::ParseUInt64(tokens[9]);

            this->dataHeader.BoundingBox().Set(minX, minY, minZ, maxX, maxY, maxZ);
            this->dataHeader.SetHasIDs(hasIDs);
            this->dataHeader.SetParticleCount(partCount);
            this->dataHeader.SetTimeCount(timeCount);
            this->dataHeader.Types().SetCount(typeCount);

        } catch(...) {
            vislib::StringA msg;
            msg.Format("Failed to parse header line file \"%s\"", fieldName);
            _ERROR_OUT(fieldName);
        }

        // reading particle types
        // Note: This is the only place where 'unicode' is really relevant!
        for (UINT32 typeIdx = 0; typeIdx < this->dataHeader.Types().Count(); typeIdx++) {
            MMSPDHeader::TypeDefinition &type = this->dataHeader.Types()[typeIdx];
            UINT32 constFieldCnt, fieldCnt;
            vislib::StringA str;

            _ASSERT_READLINE(line);
            if (unicode) {
                vislib::StringW uniLine;
                if (!vislib::UTF8Encoder::Decode(uniLine, line)) {
                    vislib::StringA msg;
                    msg.Format("Failed to decode UTF8 particle line %d", static_cast<int>(typeIdx));
                    _ERROR_OUT(msg);
                }
                line = uniLine;
            }
            line.TrimSpaces();
            tokens = vislib::StringTokeniserA::Split(line, ' ', true);
            if (tokens.Count() < 3) {
                vislib::StringA msg;
                msg.Format("Particle line %d incomplete", static_cast<int>(typeIdx));
                _ERROR_OUT(msg);
            }

            try {
                constFieldCnt = static_cast<UINT32>(vislib::CharTraitsA::ParseUInt64(tokens[1]));
            } catch(...) {
                vislib::StringA msg;
                msg.Format("Failed to parse fixFieldCount of particle line %d", static_cast<int>(typeIdx));
                _ERROR_OUT(msg);
            }

            try {
                fieldCnt = static_cast<UINT32>(vislib::CharTraitsA::ParseUInt64(tokens[2]));
            } catch(...) {
                vislib::StringA msg;
                msg.Format("Failed to parse varFieldCount of particle line %d", static_cast<int>(typeIdx));
                _ERROR_OUT(msg);
            }

            if (tokens.Count() < 3 + constFieldCnt * 3 + fieldCnt * 2) {
                vislib::StringA msg;
                msg.Format("Particle line %d incomplete", static_cast<int>(typeIdx));
                _ERROR_OUT(msg);
            } else if (tokens.Count() > 3 + constFieldCnt * 3 + fieldCnt * 2) {
                vislib::sys::Log::DefaultLog.WriteWarn("Trailing information on particle line %d will be ignored", static_cast<int>(typeIdx));
            }

            type.SetBaseType(tokens[0]);
            type.ConstFields().SetCount(constFieldCnt);
            type.Fields().SetCount(fieldCnt);
            int pos = 3;
            
            for (UINT32 fieldIdx = 0; fieldIdx < constFieldCnt; fieldIdx++, pos += 3) {
                MMSPDHeader::ConstField &field = type.ConstFields()[fieldIdx];

                if (tokens[pos].Equals("id")) _ERROR_OUT("Field \"id\" is reserved for internal use and must not be used!");
                if (tokens[pos].Equals("type")) _ERROR_OUT("Field \"type\" is reserved for internal use and must not be used!");
                if (tokens[pos + 1].Equals("b") || tokens[pos + 1].Equals("byte", false)) field.SetType(MMSPDHeader::Field::TYPE_BYTE);
                else if (tokens[pos + 1].Equals("f") || tokens[pos + 1].Equals("float", false)) field.SetType(MMSPDHeader::Field::TYPE_FLOAT);
                else if (tokens[pos + 1].Equals("d") || tokens[pos + 1].Equals("double", false)) field.SetType(MMSPDHeader::Field::TYPE_DOUBLE);
                else {
                    str.Format("Type \"%s\" of field \"%d\" of type definition \"%d\" is unknown",
                        tokens[pos + 1].PeekBuffer(), static_cast<int>(fieldIdx), static_cast<int>(typeIdx));
                   _ERROR_OUT(str);
                }
                field.SetName(tokens[pos]);
                try {
                    switch (field.GetType()) {
                        case MMSPDHeader::Field::TYPE_BYTE: {
                            int i = vislib::CharTraitsA::ParseInt(tokens[pos + 2]);
                            if ((i < 0) || (i > 255)) {
                                str.Format("Byte value of field \"%s\" of type %d out of range\n",
                                    tokens[pos].PeekBuffer(), static_cast<int>(typeIdx));
                               _ERROR_OUT(str);
                            }
                            field.SetByte(static_cast<unsigned char>(i));
                        } break;
                        case MMSPDHeader::Field::TYPE_FLOAT: {
                            field.SetFloat(static_cast<float>(vislib::CharTraitsA::ParseDouble(tokens[pos + 2])));
                        } break;
                        case MMSPDHeader::Field::TYPE_DOUBLE: {
                            field.SetDouble(vislib::CharTraitsA::ParseDouble(tokens[pos + 2]));
                        } break;
                        default: _ERROR_OUT("Internal Error!");
                    }
                } catch(...) {
                    str.Format("Failed to parse value for field \"%s\" of type %d\n",
                        tokens[pos].PeekBuffer(), static_cast<int>(typeIdx));
                   _ERROR_OUT(str);
                }
            }

            for (UINT32 fieldIdx = 0; fieldIdx < fieldCnt; fieldIdx++, pos += 2) {
                MMSPDHeader::Field &field = type.Fields()[fieldIdx];

                if (tokens[pos].Equals("id")) _ERROR_OUT("Field \"id\" is reserved for internal use and must not be used!");
                if (tokens[pos].Equals("type")) _ERROR_OUT("Field \"type\" is reserved for internal use and must not be used!");
                if (tokens[pos + 1].Equals("b") || tokens[pos + 1].Equals("byte", false)) field.SetType(MMSPDHeader::Field::TYPE_BYTE);
                else if (tokens[pos + 1].Equals("f") || tokens[pos + 1].Equals("float", false)) field.SetType(MMSPDHeader::Field::TYPE_FLOAT);
                else if (tokens[pos + 1].Equals("d") || tokens[pos + 1].Equals("double", false)) field.SetType(MMSPDHeader::Field::TYPE_DOUBLE);
                else {
                    str.Format("Type \"%s\" of field \"%d\" of type definition \"%d\" is unknown",
                        tokens[pos + 1].PeekBuffer(), static_cast<int>(fieldIdx), static_cast<int>(typeIdx));
                   _ERROR_OUT(str);
                }
                field.SetName(str);
            }

        }

    } else {
        // reading header line
        unsigned char hasIDs;
        double minX, minY, minZ, maxX, maxY, maxZ;
        UINT32 timeCount;
        UINT32 typeCount;
        UINT64 partCount;

        _ASSERT_READFILE(&hasIDs, 1);
        _ASSERT_READFILE(&minX, 8);
        _ASSERT_READFILE(&minY, 8);
        _ASSERT_READFILE(&minZ, 8);
        _ASSERT_READFILE(&maxX, 8);
        _ASSERT_READFILE(&maxY, 8);
        _ASSERT_READFILE(&maxZ, 8);
        _ASSERT_READFILE(&timeCount, 4);
        _ASSERT_READFILE(&typeCount, 4);
        _ASSERT_READFILE(&partCount, 8);

        // now I am confident enough to start setting data
        this->dataHeader.BoundingBox().Set(minX, minY, minZ, maxX, maxY, maxZ);
        this->dataHeader.SetHasIDs(hasIDs != 0);
        this->dataHeader.SetParticleCount(partCount);
        this->dataHeader.SetTimeCount(timeCount);
        this->dataHeader.Types().SetCount(typeCount);

        // reading particle types
        for (UINT32 typeIdx = 0; typeIdx < typeCount; typeIdx++) {
            MMSPDHeader::TypeDefinition &type = this->dataHeader.Types()[typeIdx];
            vislib::StringA str;
            UINT32 constFieldCnt, fieldCnt;

            _ASSERT_READSTRINGBINARY(str);
            _ASSERT_READFILE(&constFieldCnt, 4);
            if (bigEndian) {
                unsigned char *fac = reinterpret_cast<unsigned char*>(&constFieldCnt);
                vislib::Swap(fac[0], fac[3]);
                vislib::Swap(fac[1], fac[2]);
            }
            _ASSERT_READFILE(&fieldCnt, 4);
            if (bigEndian) {
                unsigned char *fac = reinterpret_cast<unsigned char*>(&fieldCnt);
                vislib::Swap(fac[0], fac[3]);
                vislib::Swap(fac[1], fac[2]);
            }

            type.SetBaseType(str);
            type.ConstFields().SetCount(constFieldCnt);
            type.Fields().SetCount(fieldCnt);

            for (UINT32 fieldIdx = 0; fieldIdx < constFieldCnt; fieldIdx++) {
                vislib::StringA typeStr;
                MMSPDHeader::ConstField &field = type.ConstFields()[fieldIdx];

                _ASSERT_READSTRINGBINARY(str);
                if (str.Equals("id")) _ERROR_OUT("Field \"id\" is reserved for internal use and must not be used!");
                if (str.Equals("type")) _ERROR_OUT("Field \"type\" is reserved for internal use and must not be used!");
                _ASSERT_READSTRINGBINARY(typeStr);
                if (typeStr.Equals("b") || typeStr.Equals("byte", false)) field.SetType(MMSPDHeader::Field::TYPE_BYTE);
                else if (typeStr.Equals("f") || typeStr.Equals("float", false)) field.SetType(MMSPDHeader::Field::TYPE_FLOAT);
                else if (typeStr.Equals("d") || typeStr.Equals("double", false)) field.SetType(MMSPDHeader::Field::TYPE_DOUBLE);
                else {
                    str.Format("Type \"%s\" of field \"%d\" of type definition \"%d\" is unknown",
                        typeStr.PeekBuffer(), static_cast<int>(fieldIdx), static_cast<int>(typeIdx));
                   _ERROR_OUT(str);
                }
                field.SetName(str);
                switch (field.GetType()) {
                    case MMSPDHeader::Field::TYPE_BYTE: {
                        unsigned char b;
                        _ASSERT_READFILE(&b, 1);
                        field.SetByte(b);
                    } break;
                    case MMSPDHeader::Field::TYPE_FLOAT: {
                        float f;
                        _ASSERT_READFILE(&f, 4);
                        if (bigEndian) {
                            unsigned char *fac = reinterpret_cast<unsigned char*>(&f);
                            vislib::Swap(fac[0], fac[3]);
                            vislib::Swap(fac[1], fac[2]);
                        }
                        field.SetFloat(f);
                    } break;
                    case MMSPDHeader::Field::TYPE_DOUBLE: {
                        double d;
                        _ASSERT_READFILE(&d, 8);
                        if (bigEndian) {
                            unsigned char *fac = reinterpret_cast<unsigned char*>(&d);
                            vislib::Swap(fac[0], fac[7]);
                            vislib::Swap(fac[1], fac[6]);
                            vislib::Swap(fac[2], fac[5]);
                            vislib::Swap(fac[3], fac[4]);
                        }
                        field.SetDouble(d);
                    } break;
                    default: _ERROR_OUT("Internal Error!");
                }
            }

            for (UINT32 fieldIdx = 0; fieldIdx < fieldCnt; fieldIdx++) {
                vislib::StringA typeStr;
                MMSPDHeader::Field &field = type.Fields()[fieldIdx];

                _ASSERT_READSTRINGBINARY(str);
                if (str.Equals("id")) _ERROR_OUT("Field \"id\" is reserved for internal use and must not be used!");
                if (str.Equals("type")) _ERROR_OUT("Field \"type\" is reserved for internal use and must not be used!");
                _ASSERT_READSTRINGBINARY(typeStr);
                if (typeStr.Equals("b") || typeStr.Equals("byte", false)) field.SetType(MMSPDHeader::Field::TYPE_BYTE);
                else if (typeStr.Equals("f") || typeStr.Equals("float", false)) field.SetType(MMSPDHeader::Field::TYPE_FLOAT);
                else if (typeStr.Equals("d") || typeStr.Equals("double", false)) field.SetType(MMSPDHeader::Field::TYPE_DOUBLE);
                else {
                    str.Format("Type \"%s\" of field \"%d\" of type definition \"%d\" is unknown",
                        typeStr.PeekBuffer(), static_cast<int>(fieldIdx), static_cast<int>(typeIdx));
                   _ERROR_OUT(str);
                }
                field.SetName(str);
            }

        }

    }

    // reading frames
    //  index generation and size estimation
    this->frameIdx = new UINT64[this->dataHeader.GetTimeCount() + 1];
    this->framePartCnts = new UINT64[this->dataHeader.GetTimeCount()];
    this->frameIdx[0] = static_cast<UINT64>(this->file->Tell());

    // TODO: Implement

/*

    UINT32 frmCnt = 0;
    _ASSERT_READFILE(&frmCnt, 4);
    if (frmCnt == 0) {
        _ERROR_OUT("MMPLD file does not contain any frame information");
    }

    float box[6];
    _ASSERT_READFILE(box, 4 * 6);
    this->bbox.Set(box[0], box[1], box[2], box[3], box[4], box[5]);
    _ASSERT_READFILE(box, 4 * 6);
    this->clipbox.Set(box[0], box[1], box[2], box[3], box[4], box[5]);

    delete[] this->frameIdx;
    this->frameIdx = new UINT64[frmCnt + 1];
    _ASSERT_READFILE(this->frameIdx, 8 * (frmCnt + 1));
    double size = 0.0;
    for (UINT32 i = 0; i < frmCnt; i++) {
        size += static_cast<double>(this->frameIdx[i + 1] - this->frameIdx[i]);
    }
    size /= static_cast<double>(frmCnt);
    size *= CACHE_FRAME_FACTOR;

    UINT64 mem = vislib::sys::SystemInformation::AvailableMemorySize();
    unsigned int cacheSize = static_cast<unsigned int>(mem / size);

    if (cacheSize > CACHE_SIZE_MAX) {
        cacheSize = CACHE_SIZE_MAX;
    }
    if (cacheSize < CACHE_SIZE_MIN) {
        vislib::StringA msg;
        msg.Format("Frame cache size forced to %i. Calculated size was %u.\n",
            CACHE_SIZE_MIN, cacheSize);
        this->GetCoreInstance()->Log().WriteMsg(vislib::sys::Log::LEVEL_WARN, msg);
        cacheSize = CACHE_SIZE_MIN;
    } else {
        vislib::StringA msg;
        msg.Format("Frame cache size set to %i.\n", cacheSize);
        this->GetCoreInstance()->Log().WriteMsg(vislib::sys::Log::LEVEL_INFO, msg);
    }

    this->setFrameCount(frmCnt);
    this->initFrameCache(cacheSize);
    */

#undef _ASSERT_READLINE
#undef _ASSERT_READSTRINGBINARY
#undef _ASSERT_READFILE
#undef _ERROR_OUT

    return true;
}


/*
 * moldyn::MMSPDDataSource::getDataCallback
 */
bool moldyn::MMSPDDataSource::getDataCallback(Call& caller) {
    MultiParticleDataCall *c2 = dynamic_cast<MultiParticleDataCall*>(&caller);
    /*if (c2 == NULL)*/ return false;

    //Frame *f = NULL;
    //if (c2 != NULL) {
    //    f = dynamic_cast<Frame *>(this->requestLockedFrame(c2->FrameID()));
    //    if (f == NULL) return false;
    //    c2->SetUnlocker(new Unlocker(*f));
    //    c2->SetFrameID(f->FrameNumber());
    //    c2->SetDataHash(0);
    //    f->SetData(*c2);
    //}

    //return true;
}


/*
 * moldyn::MMSPDDataSource::getExtentCallback
 */
bool moldyn::MMSPDDataSource::getExtentCallback(Call& caller) {
    MultiParticleDataCall *c2 = dynamic_cast<MultiParticleDataCall*>(&caller);

    if (c2 != NULL) {
        c2->SetFrameCount(vislib::math::Max(1u, this->dataHeader.GetTimeCount()));
        c2->AccessBoundingBoxes().Clear();
        c2->AccessBoundingBoxes().SetObjectSpaceBBox(this->dataHeader.GetBoundingBox());
        c2->AccessBoundingBoxes().SetObjectSpaceClipBox(this->clipbox);
        c2->SetUnlocker(NULL);
        return true;
    }

    return false;
}
