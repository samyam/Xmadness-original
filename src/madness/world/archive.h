/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

#ifndef MADNESS_WORLD_ARCHIVE_H__INCLUDED
#define MADNESS_WORLD_ARCHIVE_H__INCLUDED

/**
 \file archive.h
 \brief Interface templates for the archives (serialization).
 \ingroup serialization
*/

#include <type_traits>
#include <complex>
#include <iostream>
#include <cstdio>
#include <vector>
#include <map>
#include <tuple>
//#include <madness/world/worldprofile.h>
#include <madness/world/type_traits.h>
#include <madness/world/madness_exception.h>

/// \todo Brief description needed.
#define ARCHIVE_COOKIE "archive"

/// Major version number for archive.
#define ARCHIVE_MAJOR_VERSION 0
/// Minor version number for archive.
#define ARCHIVE_MINOR_VERSION 1

//#define MAD_ARCHIVE_DEBUG_ENABLE

/// Macro for helping debug archive tools.
#ifdef MAD_ARCHIVE_DEBUG_ENABLE
#define MAD_ARCHIVE_DEBUG(s) s
//using std::endl;
#else
#define MAD_ARCHIVE_DEBUG(s)
#endif

namespace madness {

    // Forward declarations
    template <typename T> class Tensor;

    namespace archive {

        // Forward declarations
        template <class>
        class archive_array;
        template <class T>
        inline archive_array<T> wrap(const T*, unsigned int);
        template <class T>
        inline archive_array<unsigned char> wrap_opaque(const T*, unsigned int);
        template <class T>
        inline archive_array<unsigned char> wrap_opaque(const T&);

        /// \addtogroup serialization
        /// @{

        // There are 64 empty slots for user types. Free space for
        // registering user types begins at cookie=128.



        /// The list of type names for use in archives.

        /// \todo Could this namespace-scope variable be defined in archive_type_names.cc and just externed here always?
#ifdef MAD_ARCHIVE_TYPE_NAMES_CC
        const char *archive_type_names[256];
#else
        extern const char *archive_type_names[256];
#endif

        /// Initializes the type names for the archives.
        // Implemented in archive_type_names.cc
        void archive_initialize_type_names();

        /// Used to enable type checking inside archives.

        /// \tparam T The data type.
        template <typename T>
        struct archive_typeinfo {
            static const unsigned char cookie = 255; ///< Numeric ID for the type; 255 indicates unknown type.
        };

        /// Returns the name of the type, or unknown if not registered.

        /// \tparam T The data type.
        /// \return The name of the type.
        template <typename T>
        const char* get_type_name() {
            return archive_type_names[archive_typeinfo<T>::cookie];
        }

        /// \todo Brief description needed (ARCHIVE_REGISTER_TYPE_XLC_EXTRA).

        /// \param[in] T The type to register.
#if defined(ARCHIVE_REGISTER_TYPE_INSTANTIATE_HERE) && defined(ARCHIVE_REGISTER_TYPE_IBMBUG)
#define ARCHIVE_REGISTER_TYPE_XLC_EXTRA(T) \
        ; const unsigned char archive_typeinfo< T >::cookie
#else
#define ARCHIVE_REGISTER_TYPE_XLC_EXTRA(T)
#endif

        /// Used to associate a type with a cookie value inside archive.

        /// Makes a specialization of \c archive_typeinfo for type \c T that
        /// specifies the correct cookie value.
        /// \param[in] T The type.
        /// \param[in] cooky The cookie value.
#define ARCHIVE_REGISTER_TYPE(T, cooky) \
        template <> \
        struct archive_typeinfo< T > { \
            static const unsigned char cookie = cooky; \
        } \
        ARCHIVE_REGISTER_TYPE_XLC_EXTRA(T)


        /// Used to associate a type and a pointer to the type with a cookie value inside archive.

        /// \param[in] T The type.
        /// \param[in] cooky The cookie value.
#define ARCHIVE_REGISTER_TYPE_AND_PTR(T, cooky) \
        ARCHIVE_REGISTER_TYPE(T, cooky); \
        ARCHIVE_REGISTER_TYPE(T*, cooky+64)


        /// Alias for \c archive_type_names.
#define ATN ::madness::archive::archive_type_names
        /// Alias for \c archive_typeinfo.
#define ATI ::madness::archive::archive_typeinfo


        /// Used to associate names with types.

        /// \param[in] T The type.
#define ARCHIVE_REGISTER_TYPE_NAME(T) \
        if (strcmp( ATN[ATI< T >::cookie], "invalid") ) { \
            std::cout << "archive_register_type_name: slot/cookie already in use! " << #T << " " << ATN[ATI< T >::cookie] << std::endl; \
            MADNESS_EXCEPTION("archive_register_type_name: slot/cookie already in use!", 0); \
         } \
         ATN[ATI< T >::cookie] = #T        


        /// Used to associate names with types and pointers to that type.

        /// \param[in] T The type.
#define ARCHIVE_REGISTER_TYPE_AND_PTR_NAMES(T) \
        ARCHIVE_REGISTER_TYPE_NAME(T); \
        ARCHIVE_REGISTER_TYPE_NAME(T*)



#ifndef DOXYGEN_SHOULD_SKIP_THIS
        // **********
        // Register standard types and "common" MADNESS types.
        //
        // doxygen interprets these all as functions, not as calls to a macro.
        // Thus, we force doxygen to skip this block.

        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned char,0);
        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned short,1);
        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned int,2);
        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned long,3);
        ARCHIVE_REGISTER_TYPE_AND_PTR(unsigned long long,4);
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed char,5);
        ARCHIVE_REGISTER_TYPE_AND_PTR(char,5);	// Needed, but why?
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed short,6);
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed int,7);
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed long,8);
        ARCHIVE_REGISTER_TYPE_AND_PTR(signed long long,9);
        ARCHIVE_REGISTER_TYPE_AND_PTR(bool,10);
        ARCHIVE_REGISTER_TYPE_AND_PTR(float,11);
        ARCHIVE_REGISTER_TYPE_AND_PTR(double,12);
        ARCHIVE_REGISTER_TYPE_AND_PTR(long double,13);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::complex<float>,14);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::complex<double>,15);

        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<char>,20);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<unsigned char>,21);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<short>,22);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<unsigned short>,23);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<int>,24);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<unsigned int>,25);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<long>,26);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<unsigned long>,27);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<bool>,28);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<float>,29);
        ARCHIVE_REGISTER_TYPE_AND_PTR(std::vector<double>,30);

        ARCHIVE_REGISTER_TYPE_AND_PTR(std::string,31);

        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor<int>,32);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor<long>,33);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor<float>,34);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor<double>,35);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor< std::complex<float> >,36);
        ARCHIVE_REGISTER_TYPE_AND_PTR(Tensor< std::complex<double> >,37);
        // **********
#endif


        /// Base class for all archive classes.
        class BaseArchive {
        public:
            static const bool is_archive = true; ///< Flag to determine if this object is an archive.
            static const bool is_input_archive = false; ///< Flag to determine if this object is an input archive.
            static const bool is_output_archive = false; ///< Flag to determine if this object is an output archive.
            static const bool is_parallel_archive = false; ///< Flag to determine if this object is a parallel archive.

            BaseArchive() {
                archive_initialize_type_names();
            }
        }; // class BaseArchive


        /// Base class for input archive classes.
        class BaseInputArchive : public BaseArchive {
        public:
            static const bool is_input_archive = true; ///< Flag to determine if this object is an input archive.
        }; // class BaseInputArchive


        /// Base class for output archive classes.
        class BaseOutputArchive : public BaseArchive {
        public:
            static const bool is_output_archive = true; ///< Flag to determine if this object is an output archive.
        }; // class BaseOutputArchive


        /// Checks if \c T is an archive type.

        /// If \c T is an archive type, then \c is_archive will be inherited
        /// from \c std::true_type, otherwise it is inherited from
        /// \c std::false_type.
        /// \tparam T The type to check.
        template <typename T>
        struct is_archive : public std::is_base_of<BaseArchive, T>{};


        /// Checks if \c T is an input archive type.

        /// If \c T is an input archive type, then \c is_input_archive will be
        /// inherited from \c std::true_type, otherwise it is inherited from
        /// \c std::false_type.
        /// \tparam T The type to check.
        template <typename T>
        struct is_input_archive : public std::is_base_of<BaseInputArchive, T> {};


        /// Checks if \c T is an output archive type.

        /// If \c T is an output archive type, then \c is_output_archive will
        /// be inherited from \c std::true_type, otherwise it is inherited from
        /// \c std::false_type.
        /// \tparam T The type to check.
        template <typename T>
        struct is_output_archive : public std::is_base_of<BaseOutputArchive, T> {};


        /// Serialize an array of fundamental stuff.

        /// The function only appears (via \c enable_if) if \c T is
        /// serializable and \c Archive is an output archive.
        /// \tparam Archive The archive type.
        /// \tparam T The type of data in the array.
        /// \param[in] ar The archive.
        /// \param[in] t Pointer to the start of the array.
        /// \param[in] n Number of data items to be serialized.
        template <class Archive, class T>
        typename std::enable_if< is_serializable<T>::value && is_output_archive<Archive>::value >::type
        serialize(const Archive& ar, const T* t, unsigned int n) {
            MAD_ARCHIVE_DEBUG(std::cout << "serialize fund array" << std::endl);
            ar.store(t,n);
        }


        /// Deserialize an array of fundamental stuff.

        /// The function only appears (via \c enable_if) if \c T is
        /// serializable and \c Archive is an input archive.
        /// \tparam Archive The archive type.
        /// \tparam T The type of data in the array.
        /// \param[in] ar The archive.
        /// \param[in] t Pointer to the start of the array.
        /// \param[in] n Number of data items to be deserialized.
        template <class Archive, class T>
        typename std::enable_if< is_serializable<T>::value && is_input_archive<Archive>::value >::type
        serialize(const Archive& ar, const T* t, unsigned int n) {
            MAD_ARCHIVE_DEBUG(std::cout << "deserialize fund array" << std::endl);
            ar.load((T*) t,n);
        }


        /// Serialize (or deserialize) an array of non-fundamental stuff.

        /// The function only appears (via \c enable_if) if \c T is
        /// not serializable and \c Archive is an archive.
        /// \tparam Archive The archive type.
        /// \tparam T The type of data in the array.
        /// \param[in] ar The archive.
        /// \param[in] t Pointer to the start of the array.
        /// \param[in] n Number of data items to be serialized.
        template <class Archive, class T>
        typename std::enable_if< ! is_serializable<T>::value && is_archive<Archive>::value >::type
        serialize(const Archive& ar, const T* t, unsigned int n) {
            MAD_ARCHIVE_DEBUG(std::cout << "(de)serialize non-fund array" << std::endl);
            for (unsigned int i=0; i<n; ++i)
                ar & t[i];
        }


        /// Default implementation of the pre/postamble for type checking.

        /// \tparam Archive The archive class.
        /// \tparam T The type to serialized or to expect upon deserialization.
        template <class Archive, class T>
        struct ArchivePrePostImpl {
            /// Deserialize a cookie and check the type.

            /// \param[in] ar The archive.
            static inline void preamble_load(const Archive& ar) {
                unsigned char ck = archive_typeinfo<T>::cookie;
                unsigned char cookie;
                ar.load(&cookie, 1); // cannot use >>
                if (cookie != ck) {
                    char msg[255];
                    std::sprintf(msg,"InputArchive type mismatch: expected cookie "
                                 "%u (%s) but got %u (%s) instead",
                                 ck, archive_type_names[ck],
                                 cookie,archive_type_names[cookie]);
                    std::cerr << msg << std::endl;
                    MADNESS_EXCEPTION(msg, static_cast<int>(cookie));
                }
                else {
                    MAD_ARCHIVE_DEBUG(std::cout << "read cookie " << archive_type_names[cookie] << std::endl);
                }
            }

            /// Serialize a cookie for type checking.

            /// \param[in] ar The archive.
            static inline void preamble_store(const Archive& ar) {
                unsigned char ck = archive_typeinfo<T>::cookie;
                ar.store(&ck, 1); // cannot use <<
                MAD_ARCHIVE_DEBUG(std::cout << "wrote cookie " << archive_type_names[ck] << std::endl);
            }

            /// By default there is no postamble.
            static inline void postamble_load(const Archive& /*ar*/) {}

            /// By default there is no postamble.
            static inline void postamble_store(const Archive& /*ar*/) {}
        };


        /// Default symmetric serialization of a non-fundamental type.

        /// \tparam Archive The archive type.
        /// \tparam T The type to symmetrically serialize.
        template <class Archive, class T>
        struct ArchiveSerializeImpl {
            /// Serializes the type.

            /// \param[in] ar The archive.
            /// \param[in,out] t The data.
            static inline void serialize(const Archive& ar, T& t) {
                t.serialize(ar);
            }
        };


        /// Redirect `serialize(ar, t)` to `serialize(ar, &t, 1)` for fundamental types.

        /// The function only appears (due to \c enable_if) if \c T is
        /// serializable and \c Archive is an archive.
        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        /// \param[in] ar The archive.
        /// \param[in] t The data to be serialized.
        template <class Archive, class T>
        inline
        typename std::enable_if< is_serializable<T>::value && is_archive<Archive>::value >::type
        serialize(const Archive& ar, const T& t) {
            MAD_ARCHIVE_DEBUG(std::cout << "serialize(ar,t) -> serialize(ar,&t,1)" << std::endl);
            serialize(ar,&t,1);
        }


        /// Redirect `serialize(ar,t)` to \c ArchiveSerializeImpl for non-fundamental types.

        /// The function only appears (due to \c enable_if) if \c T is not
        /// serializable and \c Archive is an archive.
        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        /// \param[in] ar The archive.
        /// \param[in] t The data to be serialized.
        template <class Archive, class T>
        inline
        typename std::enable_if< !is_serializable<T>::value && is_archive<Archive>::value >::type
        serialize(const Archive& ar, const T& t) {
            MAD_ARCHIVE_DEBUG(std::cout << "serialize(ar,t) -> ArchiveSerializeImpl" << std::endl);
            ArchiveSerializeImpl<Archive,T>::serialize(ar,(T&) t);
        }


        /// Default store of an object via `serialize(ar, t)`.

        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        template <class Archive, class T>
        struct ArchiveStoreImpl {
            /// Store an object.

            /// \param[in] ar The archive.
            /// \param[in] t The data.
            static inline void store(const Archive& ar, const T& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "store(ar,t) default" << std::endl);
                serialize(ar,t);
            }
        };


        /// Default load of an object via `serialize(ar, t)`.

        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        template <class Archive, class T>
        struct ArchiveLoadImpl {
            /// Load an object.

            /// \param[in] ar The archive.
            /// \param[in] t The data.
            static inline void load(const Archive& ar, const T& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "load(ar,t) default" << std::endl);
                serialize(ar,t);
            }
        };


        /// Default implementations of \c wrap_store and \c wrap_load.

        /// "Wrapping" refers to the addition of the type's preamble and
        /// postamble around the data to provide runtime type-checking.
        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        template <class Archive, class T>
        struct ArchiveImpl {
            /// Store an object sandwiched between its preamble and postamble.

            /// \param[in] ar The archive.
            /// \param[in] t The data.
            /// \return The archive.
            static inline const Archive& wrap_store(const Archive& ar, const T& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_store for default" << std::endl);
                ArchivePrePostImpl<Archive,T>::preamble_store(ar);
                ArchiveStoreImpl<Archive,T>::store(ar,t);
                ArchivePrePostImpl<Archive,T>::postamble_store(ar);
                return ar;
            }

            /// Load an object sandwiched between its preamble and postamble.

            /// \param[in] ar The archive.
            /// \param[in] t The data.
            /// \return The archive.
            static inline const Archive& wrap_load(const Archive& ar, const T& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_load for default" << std::endl);
                ArchivePrePostImpl<Archive,T>::preamble_load(ar);
                ArchiveLoadImpl<Archive,T>::load(ar,(T&) t);  // Loses constness here!
                ArchivePrePostImpl<Archive,T>::postamble_load(ar);
                return ar;
            }
        };


        /// Redirect \c << to \c ArchiveImpl::wrap_store for output archives.

        /// The function only appears (due to \c enable_if) if \c Archive
        /// is an output archive.
        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        /// \param[in] ar The archive.
        /// \param[in] t The data.
        template <class Archive, class T>
        inline
        typename std::enable_if<is_output_archive<Archive>::value, const Archive&>::type
        operator<<(const Archive& ar, const T& t) {
            //PROFILE_FUNC;
            return ArchiveImpl<Archive,T>::wrap_store(ar,t);
        }

        /// Redirect \c >> to `ArchiveImpl::wrap_load` for input archives.

        /// The function only appears (due to \c enable_if) if \c Archive
        /// is an input archive.
        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        /// \param[in] ar The archive.
        /// \param[in] t The data.
        template <class Archive, class T>
        inline
        typename std::enable_if<is_input_archive<Archive>::value, const Archive&>::type
        operator>>(const Archive& ar, const T& t) {
            //PROFILE_FUNC;
            return ArchiveImpl<Archive,T>::wrap_load(ar,t);
        }

        /// Redirect \c & to `ArchiveImpl::wrap_store` for output archives.

        /// The function only appears (due to \c enable_if) if \c Archive
        /// is an output archive.
        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        /// \param[in] ar The archive.
        /// \param[in] t The data.
        template <class Archive, class T>
        inline
        typename std::enable_if<is_output_archive<Archive>::value, const Archive&>::type
        operator&(const Archive& ar, const T& t) {
            //PROFILE_FUNC;
            return ArchiveImpl<Archive,T>::wrap_store(ar,t);
        }

        /// Redirect \c & to `ArchiveImpl::wrap_load` for input archives.

        /// The function only appears (due to \c enable_if) if \c Archive
        /// is an input archive.
        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        /// \param[in] ar The archive.
        /// \param[in] t The data.
        template <class Archive, class T>
        inline
        typename std::enable_if<is_input_archive<Archive>::value, const Archive&>::type
        operator&(const Archive& ar, const T& t) {
            //PROFILE_FUNC;
            return ArchiveImpl<Archive,T>::wrap_load(ar,t);
        }


        // -----------------------------------------------------------------

        /// Wrapper for an opaque pointer for serialization purposes.

        /// Performs a bitwise copy of the pointer without any remapping.
        /// \tparam T The type of object being pointed to.
        /// \todo Verify this documentation.
        template <class T>
        class archive_ptr {
        public:
            T* ptr; ///< The pointer.

            /// Constructor specifying \c nullptr by default.

            /// \param[in] t The pointer.
            archive_ptr(T* t = nullptr)
                : ptr(t) {}

            /// Dereference the pointer.

            /// \return The dereferenced pointer.
            T& operator*() {
                return *ptr;
            }

            /// Serialize the pointer.

            /// \tparam Archive The archive type.
            /// \param[in] ar The archive.
            template <class Archive>
            void serialize(const Archive& ar) {ar & wrap_opaque(&ptr, 1);}
        };


        /// Wrapper for pointers.

        /// \tparam T The type of object being pointed to.
        /// \param[in] p The pointer.
        /// \return The wrapped pointer.
      	template <class T>
      	inline archive_ptr<T> wrap_ptr(T* p) {
            return archive_ptr<T>(p);
        }

        /// Wrapper for dynamic arrays and pointers.

        /// \tparam T The type of object being pointed to.
        template <class T>
        class archive_array {
        public:
            const T* ptr; ///< The pointer.
            unsigned int n; ///< The number of objects in the array.

            /// Constructor specifying a memory location and size.

            /// \param[in] ptr The pointer.
            /// \param[in] n The number of objects in the array.
            archive_array(const T *ptr, unsigned int n) : ptr(ptr), n(n) {}

            /// Constructor specifying no array and of 0 length.
            archive_array() : ptr(nullptr), n(0) {}
        };


        /// Factory function to wrap a dynamically allocated pointer as a typed \c archive_array.

        /// \tparam T The data type.
        /// \param[in] ptr The pointer.
        /// \param[in] n The number of data elements in the array.
        /// \return The wrapped pointer.
        template <class T>
        inline archive_array<T> wrap(const T* ptr, unsigned int n) {
            return archive_array<T>(ptr,n);
        }


        /// Factory function to wrap a pointer to contiguous data as an opaque (\c uchar) \c archive_array.

        /// \tparam T The data type.
        /// \param[in] ptr The pointer.
        /// \param[in] n The number of data elements in the array.
        /// \return The wrapped pointer, as an opaque \c archive_array.
        template <class T>
        inline archive_array<unsigned char> wrap_opaque(const T* ptr, unsigned int n) {
            return archive_array<unsigned char>((unsigned char*) ptr, n*sizeof(T));
        }

        /// Factory function to wrap a contiguous scalar as an opaque (\c uchar) \c archive_array.

        /// \tparam T The data type.
        /// \param[in] t The data.
        /// \return The wrapped data.
        template <class T>
        inline archive_array<unsigned char> wrap_opaque(const T& t) {
            return archive_array<unsigned char>((unsigned char*) &t,sizeof(t));
        }


        /// Serialize a function pointer.

        /// \tparam Archive The archive type.
        /// \tparam resT The function's return type.
        /// \tparam paramT Parameter pack for the function's arguments.
        template <class Archive, typename resT, typename... paramT>
        struct ArchiveSerializeImpl<Archive, resT(*)(paramT...)> {
            /// Serialize the function pointer.

            /// \param[in] ar The archive.
            /// \param[in] fn The function pointer.
            static inline void serialize(const Archive& ar, resT(*fn)(paramT...)) {
                ar & wrap_opaque(fn);
            }
        };


        /// Serialize a member function pointer.

        /// \tparam Archive The archive type.
        /// \tparam resT The member function's return type.
        /// \tparam objT The object type.
        /// \tparam paramT Parameter pack for the member function's arguments.
        template <class Archive, typename resT, typename objT, typename... paramT>
        struct ArchiveSerializeImpl<Archive, resT(objT::*)(paramT...)> {
            /// Serialize the member function pointer.

            /// \param[in] ar The archive.
            /// \param[in] memfn The member function pointer.
            static inline void serialize(const Archive& ar, resT(objT::*memfn)(paramT...)) {
                ar & wrap_opaque(memfn);
            }
        };

        /// Serialize a const member function pointer.

        /// \tparam Archive The archive type.
        /// \tparam resT The const member function's return type.
        /// \tparam objT The object type.
        /// \tparam paramT Parameter pack for the const member function's arguments.
        template <class Archive, typename resT, typename objT, typename... paramT>
        struct ArchiveSerializeImpl<Archive, resT(objT::*)(paramT...) const> {
            /// Serialize the const member function pointer.

            /// \param[in] ar The archive.
            /// \param[in] memfn The const member function pointer.
            static inline void serialize(const Archive& ar, resT(objT::*memfn)(paramT...) const) {
                ar & wrap_opaque(memfn);
            }
        };


        /// Partial specialization of \c ArchiveImpl for \c archive_array.

        /// \tparam Archive The archive type.
        /// \tparam T The data type in the \c archive_array.
        template <class Archive, class T>
        struct ArchiveImpl< Archive, archive_array<T> > {
            /// Store the \c archive_array, wrapped by the preamble/postamble.

            /// \param[in] ar The archive.
            /// \param[in] t The \c archive_array.
            /// \return The archive.
            static inline const Archive& wrap_store(const Archive& ar, const archive_array<T>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_store for archive_array" << std::endl);
                ArchivePrePostImpl<Archive,T*>::preamble_store(ar);
                //ar << t.n;
                //ArchivePrePostImpl<Archive,T>::preamble_store(ar);
                serialize(ar,(T *) t.ptr,t.n);
                //ArchivePrePostImpl<Archive,T>::postamble_store(ar);
                ArchivePrePostImpl<Archive,T*>::postamble_store(ar);
                return ar;
            }

            /// Load the \c archive_array, using the preamble and postamble to perform runtime type-checking.

            /// \param[in] ar The archive.
            /// \param[out] t The \c archive_array.
            /// \return The archive.
            static inline const Archive& wrap_load(const Archive& ar, const archive_array<T>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_load for archive_array" << std::endl);
                ArchivePrePostImpl<Archive,T*>::preamble_load(ar);
                //unsigned int n;
                //ar >> n;
                //if (n != t.n)
                //    MADNESS_EXCEPTION("deserializing archive_array: dimension mismatch", n);
                //ArchivePrePostImpl<Archive,T>::preamble_load(ar);
                serialize(ar,(T *) t.ptr,t.n);
                //ArchivePrePostImpl<Archive,T>::postamble_load(ar);
                ArchivePrePostImpl<Archive,T*>::postamble_load(ar);
                return ar;
            }
        };


        /// Partial specialization of \c ArchiveImpl for fixed-dimension arrays that redirects to \c archive_array.

        /// \tparam Archive The archive type.
        /// \tparam T The data type.
        /// \tparam n The array size.
        template <class Archive, class T, std::size_t n>
        struct ArchiveImpl<Archive, T[n]> {
            /// Store the array, wrapped by the preamble/postamble.

            /// \param[in] ar The archive.
            /// \param[in] t The array.
            /// \return The archive.
            static inline const Archive& wrap_store(const Archive& ar, const T(&t)[n]) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_store for array" << std::endl);
                ar << wrap(&t[0],n);
                return ar;
            }

            /// Load the array, using the preamble and postamble to perform runtime type-checking.

            /// \param[in] ar The archive.
            /// \param[out] t The array.
            /// \return The archive.
            static inline const Archive& wrap_load(const Archive& ar, const T(&t)[n]) {
                MAD_ARCHIVE_DEBUG(std::cout << "wrap_load for array" << std::endl);
                ar >> wrap(&t[0],n);
                return ar;
            }
        };


        /// Serialize a complex number.

        /// \tparam Archive The archive type.
        /// \tparam T The data type underlying the complex number.
        template <class Archive, typename T>
        struct ArchiveStoreImpl< Archive, std::complex<T> > {
            /// Store a complex number.

            /// \param[in] ar The archive.
            /// \param[in] c The complex number.
            static inline void store(const Archive& ar, const std::complex<T>& c) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize complex number" << std::endl);
                ar & c.real() & c.imag();
            }
        };


        /// Deserialize a complex number.

        /// \tparam Archive the archive type.
        /// \tparam T The data type underlying the complex number.
        template <class Archive, typename T>
        struct ArchiveLoadImpl< Archive, std::complex<T> > {
            /// Load a complex number.

            /// \param[in] ar The archive.
            /// \param[out] c The complex number.
            static inline void load(const Archive& ar, std::complex<T>& c) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize complex number" << std::endl);
                T r = 0, i = 0;
                ar & r & i;
                c = std::complex<T>(r,i);
            }
        };


        /// Serialize a STL \c vector.

        /// \tparam Archive the archive type.
        /// \tparam T The data type stored in the \c vector.
        template <class Archive, typename T>
        struct ArchiveStoreImpl< Archive, std::vector<T> > {
            /// Store a \c vector.

            /// \param[in] ar The archive.
            /// \param[in] v The \c vector.
            static inline void store(const Archive& ar, const std::vector<T>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize STL vector" << std::endl);
                ar & v.size();
                ar & wrap(v.data(),v.size());
            }
        };


        /// Deserialize a STL \c vector. Clears and resizes as necessary.

        /// \tparam Archive the archive type.
        /// \tparam T The data type stored in the \c vector.
        template <class Archive, typename T>
        struct ArchiveLoadImpl< Archive, std::vector<T> > {
            /// Load a \c vector.

            /// Clears and resizes the \c vector as necessary.
            /// \param[in] ar The archive.
            /// \param[out] v The \c vector.
            static void load(const Archive& ar, std::vector<T>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize STL vector" << std::endl);
                std::size_t n = 0ul;
                ar & n;
                if (n != v.size()) {
                    v.clear();
                    v.resize(n);
                }
                ar & wrap((T *) v.data(),n);
            }
        };


        /// Serialize a STL \c vector<bool> (as a plain array of bool).

        /// \tparam Archive The archive type.
        template <class Archive>
        struct ArchiveStoreImpl< Archive, std::vector<bool> > {
            /// Store a \c vector<bool>.

            /// \param[in] ar The archive.
            /// \param[in] v The \c vector.
            static inline void store(const Archive& ar, const std::vector<bool>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize STL vector<bool>" << std::endl);
                std::size_t n = v.size();
                bool* b = new bool[n];
                for (std::size_t i=0; i<n; ++i) b[i] = v[i];
                ar & n & wrap(b,v.size());
                delete [] b;
            }
        };


        /// Deserialize a STL vector<bool>. Clears and resizes as necessary.

        /// \tparam Archive The archive type.
        template <class Archive>
        struct ArchiveLoadImpl< Archive, std::vector<bool> > {
            /// Load a \c vector<bool>.

            /// Clears and resizes the \c vector as necessary.
            /// \param[in] ar The archive.
            /// \param[out] v The \c vector.
            static void load(const Archive& ar, std::vector<bool>& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize STL vector" << std::endl);
                std::size_t n = 0ul;
                ar & n;
                if (n != v.size()) {
                    v.clear();
                    v.resize(n);
                }
                bool* b = new bool[n];
                ar & wrap(b,v.size());
                for (std::size_t i=0; i<n; ++i) v[i] = b[i];
                delete [] b;
            }
        };


        /// Serialize a STL string.

        /// \tparam Archive The archive type.
        template <class Archive>
        struct ArchiveStoreImpl< Archive, std::string > {
            /// Store a string.

            /// \param[in] ar The archive.
            /// \param[in] v The string.
            static void store(const Archive& ar, const std::string& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize STL string" << std::endl);
                ar & v.size();
                ar & wrap((const char*) v.data(),v.size());
            }
        };


        /// Deserialize a STL string. Clears and resizes as necessary.

        /// \tparam Archive The archive type.
        template <class Archive>
        struct ArchiveLoadImpl< Archive, std::string > {
            /// Load a string.

            /// Clears and resizes the string as necessary.
            /// \param[in] ar The archive.
            /// \param[out] v The string.
            static void load(const Archive& ar, std::string& v) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize STL string" << std::endl);
                std::size_t n = 0ul;
                ar & n;
                if (n != v.size()) {
                    v.clear();
                    v.resize(n);
                }
                ar & wrap((char*) v.data(),n);
            }
        };


        /// Serialize (deserialize) an STL pair.

        /// \tparam Archive The archive type.
        /// \tparam T The first data type in the pair.
        /// \tparam Q The second data type in the pair.
        template <class Archive, typename T, typename Q>
        struct ArchiveSerializeImpl< Archive, std::pair<T, Q> > {
            /// Serialize the \c pair.

            /// \param[in] ar The archive.
            /// \param[in,out] t The \c pair.
            static inline void serialize(const Archive& ar, std::pair<T,Q>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "(de)serialize STL pair" << std::endl);
                ar & t.first & t.second;
            }
        };

        namespace {

          template <size_t idx, class Archive, typename... Types>
          struct tuple_serialize_helper;

          template <class Archive, typename... Types>
          struct tuple_serialize_helper<0,Archive,Types...> {
        	  static void exec(const Archive& ar, std::tuple<Types...>& t) {
        		  ar & std::get<0>(t);
        	  }
          };
          template <size_t idx, class Archive, typename... Types>
          struct tuple_serialize_helper {
        	  static void exec(const Archive& ar, std::tuple<Types...>& t) {
        		  ar & std::get<idx>(t);
        		  tuple_serialize_helper<idx-1,Archive,Types...>::exec(ar,t);
        	  }
          };

        };

        /// Serialize (deserialize) a std::tuple

        /// \tparam Archive The archive type.
        /// \tparam Types The tuple payload
        template <class Archive, typename... Types>
        struct ArchiveSerializeImpl< Archive, std::tuple<Types...> > {
            /// Serialize the \c std::tuple.

            /// \param[in] ar The archive.
            /// \param[in,out] t The \c tuple.
            static inline void serialize(const Archive& ar, std::tuple<Types...>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "(de)serialize std::tuple" << std::endl);
                constexpr auto size = std::tuple_size<std::tuple<Types...>>::value;
                tuple_serialize_helper<size-1,Archive,Types...>::exec(ar, t);
            }
        };

        /// Serialize an STL \c map (crudely).

        /// \tparam Archive The archive type.
        /// \tparam T The map's key type.
        /// \tparam Q The map's data type.
        template <class Archive, typename T, typename Q>
        struct ArchiveStoreImpl< Archive, std::map<T,Q> > {
            /// Store a \c map.

            /// \param[in] ar The archive.
            /// \param[in] t The \c map.
            static void store(const Archive& ar, const std::map<T,Q>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "serialize STL map" << std::endl);
                ar << t.size();
                for (typename std::map<T,Q>::const_iterator p = t.begin();
                        p != t.end(); ++p) {
                    // Fun and games here since IBM's iterator (const or
                    // otherwise) gives us a const qualified key
                    // (p->first) which buggers up the type matching
                    // unless the user defines pair(T,Q) and pair(const
                    // T,Q) to have cookie (which is tedious).
                    std::pair<T,Q> pp = *p;
                    ar & pp;
                }
            }
        };


        /// Deserialize an STL \c map. The \c map is \em not cleared; duplicate elements are replaced.

        /// \tparam Archive The archive type.
        /// \tparam T The map's key type.
        /// \tparam Q The map's data type.
        template <class Archive, typename T, typename Q>
        struct ArchiveLoadImpl< Archive, std::map<T,Q> > {
            /// Load a \c map.

            /// The \c map is \em not cleared; duplicate elements are replaced.
            /// \param[in] ar The archive.
            /// \param[out] t The \c map.
            static void load(const Archive& ar, std::map<T,Q>& t) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserialize STL map" << std::endl);
                std::size_t n = 0;
                ar & n;
                while (n--) {
                    std::pair<T,Q> p;
                    ar & p;
                    t[p.first] = p.second;
                }
            }
        };

        /// @}

    }
}

#endif // MADNESS_WORLD_ARCHIVE_H__INCLUDED
