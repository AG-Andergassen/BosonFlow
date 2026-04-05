#pragma once

#include <hdf5.h>

#include <cstdio>
#include <stdexcept>
#include <string>
#include <utility>

namespace H5 {

using H5std_string = std::string;

class Exception : public std::runtime_error {
 public:
    explicit Exception(const std::string &message)
        : std::runtime_error(message) {}

    void printErrorStack() const
    {
        H5Eprint2(H5E_DEFAULT, stderr);
    }

    static void dontPrint()
    {
        H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    }
};

inline void throw_if_invalid(hid_t id, const std::string &context)
{
    if (id < 0) {
        throw Exception(context);
    }
}

inline void throw_if_error(herr_t status, const std::string &context)
{
    if (status < 0) {
        throw Exception(context);
    }
}

class Handle {
 public:
    Handle() = default;

    Handle(hid_t id, herr_t (*closer)(hid_t), bool owns = true)
        : id_(id), closer_(closer), owns_(owns) {}

    Handle(const Handle &) = delete;
    Handle &operator=(const Handle &) = delete;

    Handle(Handle &&other) noexcept
        : id_(std::exchange(other.id_, H5I_INVALID_HID)),
          closer_(std::exchange(other.closer_, nullptr)),
          owns_(std::exchange(other.owns_, false)) {}

    Handle &operator=(Handle &&other) noexcept
    {
        if (this != &other) {
            close();
            id_ = std::exchange(other.id_, H5I_INVALID_HID);
            closer_ = std::exchange(other.closer_, nullptr);
            owns_ = std::exchange(other.owns_, false);
        }
        return *this;
    }

    virtual ~Handle()
    {
        close();
    }

    hid_t raw_id() const
    {
        return id_;
    }

    void close()
    {
        if (owns_ && closer_ != nullptr && id_ >= 0) {
            closer_(id_);
        }
        id_ = H5I_INVALID_HID;
        closer_ = nullptr;
        owns_ = false;
    }

 protected:
    hid_t id_ = H5I_INVALID_HID;
    herr_t (*closer_)(hid_t) = nullptr;
    bool owns_ = false;
};

class DataType : public Handle {
 public:
    DataType() = default;

    explicit DataType(hid_t id, bool owns = false)
        : Handle(id, H5Tclose, owns) {}

    size_t getSize() const
    {
        return H5Tget_size(id_);
    }
};

struct PredType {
    inline static const DataType NATIVE_DOUBLE{H5T_NATIVE_DOUBLE, false};
    inline static const DataType NATIVE_INT{H5T_NATIVE_INT, false};
    inline static const DataType NATIVE_UINT{H5T_NATIVE_UINT, false};
    inline static const DataType C_S1{H5T_C_S1, false};
};

class StrType : public DataType {
 public:
    StrType(const DataType &base_type, size_t size)
        : DataType(H5Tcopy(base_type.raw_id()), true)
    {
        throw_if_invalid(raw_id(), "Failed to copy HDF5 string datatype");
        throw_if_error(H5Tset_size(raw_id(), size), "Failed to size HDF5 string datatype");
    }
};

class DataSpace : public Handle {
 public:
    DataSpace() = default;

    DataSpace(int rank, const hsize_t *dims)
        : Handle(H5Screate_simple(rank, dims, nullptr), H5Sclose, true)
    {
        throw_if_invalid(raw_id(), "Failed to create HDF5 dataspace");
    }

    explicit DataSpace(H5S_class_t space_type)
        : Handle(H5Screate(space_type), H5Sclose, true)
    {
        throw_if_invalid(raw_id(), "Failed to create scalar HDF5 dataspace");
    }

    explicit DataSpace(hid_t id, bool owns = true)
        : Handle(id, H5Sclose, owns) {}

    int getSimpleExtentDims(hsize_t *dims, hsize_t *maxdims) const
    {
        return H5Sget_simple_extent_dims(raw_id(), dims, maxdims);
    }
};

class DSetCreatPropList : public Handle {
 public:
    DSetCreatPropList()
        : Handle(H5Pcreate(H5P_DATASET_CREATE), H5Pclose, true)
    {
        throw_if_invalid(raw_id(), "Failed to create HDF5 dataset property list");
    }

    void setChunk(int rank, const hsize_t *dims)
    {
        throw_if_error(H5Pset_chunk(raw_id(), rank, dims), "Failed to set HDF5 chunking");
    }

    void setDeflate(unsigned level)
    {
        throw_if_error(H5Pset_deflate(raw_id(), level), "Failed to set HDF5 deflate level");
    }
};

class Attribute : public Handle {
 public:
    Attribute() = default;

    explicit Attribute(hid_t id)
        : Handle(id, H5Aclose, true) {}

    void write(const DataType &type, const void *data) const
    {
        throw_if_error(H5Awrite(raw_id(), type.raw_id(), data), "Failed to write HDF5 attribute");
    }

    void read(const DataType &type, void *data) const
    {
        throw_if_error(H5Aread(raw_id(), type.raw_id(), data), "Failed to read HDF5 attribute");
    }

    DataType getDataType() const
    {
        hid_t type_id = H5Aget_type(raw_id());
        throw_if_invalid(type_id, "Failed to get HDF5 attribute datatype");
        return DataType(type_id, true);
    }
};

class DataSet : public Handle {
 public:
    DataSet() = default;

    explicit DataSet(hid_t id)
        : Handle(id, H5Dclose, true) {}

    void write(const void *data, const DataType &type) const
    {
        throw_if_error(H5Dwrite(raw_id(), type.raw_id(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data), "Failed to write HDF5 dataset");
    }

    void read(void *data, const DataType &type) const
    {
        throw_if_error(H5Dread(raw_id(), type.raw_id(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data), "Failed to read HDF5 dataset");
    }

    DataSpace getSpace() const
    {
        hid_t space_id = H5Dget_space(raw_id());
        throw_if_invalid(space_id, "Failed to get HDF5 dataset dataspace");
        return DataSpace(space_id, true);
    }
};

class Group;

class Location : public Handle {
 public:
    using Handle::Handle;

    Group createGroup(const std::string &name) const;
    Group openGroup(const std::string &name) const;
    DataSet createDataSet(const std::string &name, const DataType &type, const DataSpace &space) const;
    DataSet createDataSet(const std::string &name, const DataType &type, const DataSpace &space, const DSetCreatPropList &plist) const;
    DataSet openDataSet(const std::string &name) const;
    Attribute createAttribute(const std::string &name, const DataType &type, const DataSpace &space) const;
    Attribute openAttribute(const std::string &name) const;
};

class Group : public Location {
 public:
    Group() = default;

    explicit Group(hid_t id)
        : Location(id, H5Gclose, true) {}
};

class H5File : public Location {
 public:
    H5File(const H5std_string &file_name, unsigned flags)
        : Location(open_file(file_name, flags), H5Fclose, true) {}

 private:
    static hid_t open_file(const H5std_string &file_name, unsigned flags)
    {
        hid_t file_id = H5I_INVALID_HID;
        if (flags == H5F_ACC_RDONLY || flags == H5F_ACC_RDWR) {
            file_id = H5Fopen(file_name.c_str(), flags, H5P_DEFAULT);
        } else {
            file_id = H5Fcreate(file_name.c_str(), flags, H5P_DEFAULT, H5P_DEFAULT);
        }
        throw_if_invalid(file_id, "Failed to open HDF5 file: " + file_name);
        return file_id;
    }
};

inline Group Location::createGroup(const std::string &name) const
{
    hid_t group_id = H5Gcreate2(raw_id(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    throw_if_invalid(group_id, "Failed to create HDF5 group: " + name);
    return Group(group_id);
}

inline Group Location::openGroup(const std::string &name) const
{
    hid_t group_id = H5Gopen2(raw_id(), name.c_str(), H5P_DEFAULT);
    throw_if_invalid(group_id, "Failed to open HDF5 group: " + name);
    return Group(group_id);
}

inline DataSet Location::createDataSet(const std::string &name, const DataType &type, const DataSpace &space) const
{
    hid_t dataset_id = H5Dcreate2(raw_id(), name.c_str(), type.raw_id(), space.raw_id(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    throw_if_invalid(dataset_id, "Failed to create HDF5 dataset: " + name);
    return DataSet(dataset_id);
}

inline DataSet Location::createDataSet(const std::string &name, const DataType &type, const DataSpace &space, const DSetCreatPropList &plist) const
{
    hid_t dataset_id = H5Dcreate2(raw_id(), name.c_str(), type.raw_id(), space.raw_id(), H5P_DEFAULT, plist.raw_id(), H5P_DEFAULT);
    throw_if_invalid(dataset_id, "Failed to create HDF5 dataset with properties: " + name);
    return DataSet(dataset_id);
}

inline DataSet Location::openDataSet(const std::string &name) const
{
    hid_t dataset_id = H5Dopen2(raw_id(), name.c_str(), H5P_DEFAULT);
    throw_if_invalid(dataset_id, "Failed to open HDF5 dataset: " + name);
    return DataSet(dataset_id);
}

inline Attribute Location::createAttribute(const std::string &name, const DataType &type, const DataSpace &space) const
{
    hid_t attribute_id = H5Acreate2(raw_id(), name.c_str(), type.raw_id(), space.raw_id(), H5P_DEFAULT, H5P_DEFAULT);
    throw_if_invalid(attribute_id, "Failed to create HDF5 attribute: " + name);
    return Attribute(attribute_id);
}

inline Attribute Location::openAttribute(const std::string &name) const
{
    hid_t attribute_id = H5Aopen(raw_id(), name.c_str(), H5P_DEFAULT);
    throw_if_invalid(attribute_id, "Failed to open HDF5 attribute: " + name);
    return Attribute(attribute_id);
}

} // namespace H5

using H5std_string = H5::H5std_string;