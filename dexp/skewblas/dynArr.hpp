// #pragma once

// #include <cstdlib>
// #include <cstring>
// #include <new>
// #include <type_traits>

// #include "arrVec.hpp"

// const INTE_TYPE _DYNAMIC_ARRAY_CONTAINER = 100;
// namespace DynamicArray
// {
//     enum class DynamicElementType
//     {
//         StasticValue,
//         DynamicValue,
//         DynamicPointer
//     };

//     template <class ELEM_TYPE>
//     class DynamicArr : public ArrVec<ELEM_TYPE>
//     {
//         // The pointer in DynamicArr is managed by the malloc~free standard.
//         // Entry-wise assignment/construction and the deconstructor are not assumed.
//         // Raw bit operations are performed in the deep copy and release actions.
//         // For usage of pointer type DynamicArr, if it owns the pointed objects,
//         // !!! manual_release_pointers() MUST BE CALLED EXPILICITLY. !!!

//         typedef DynamicArr<ELEM_TYPE> SELF_TYPE;

//     public:
//         INTE_TYPE n;
//         INTE_TYPE c;
//         DynamicElementType DynType;
//         DynamicArr() : ArrVec<ELEM_TYPE>(), n(0), c(_DYNAMIC_ARRAY_CONTAINER), DynType(StasticValue) {};
//         DynamicArr(INTE_TYPE dim, INTE_TYPE container, DynamicElementType etype) : ArrVec<ELEM_TYPE>(), n(0), c(container), DynType(etype)
//         {
//             // assert(c != 0);
//             this->d = dim;
//             if (this->d)
//             {
//                 this->v = static_cast<ELEM_TYPE *>(std::malloc(this->d * sizeof(ELEM_TYPE)));
//                 if (!this->v)
//                     throw std::bad_alloc();
//             }
//             else
//                 this->v = nullptr;
//         };
//         DynamicArr(INTE_TYPE dim, DynamicElementType etype) : DynamicArr(dim, _DYNAMIC_ARRAY_CONTAINER, etype) {};
//         void copy(const SELF_TYPE &src)
//         {
//             if (DynType != DynamicElementType::StasticValue)
//                 this->release_entry(this->n);
//             this->n = src.n;
//             this->c = src.c;
//             if (this->d == src.d)
//                 memcpy(this->v, src.v, sizeof(ELEM_TYPE) * src.n);
//             else
//             {
//                 this->release_array();
//                 this->d = src.d;
//                 this->v = static_cast<ELEM_TYPE *>(std::malloc(this->d * sizeof(ELEM_TYPE)));
//                 if (!this->v)
//                     throw std::bad_alloc();
//                 // Note that src is expected to be a legal DynamicArr that has proper nullptr initialization.
//                 // When ELEM_TYPE is a pointer type, the initalization should be carried to the new array.
//                 memcpy(this->v, src.v, sizeof(ELEM_TYPE) * src.d);
//             }
//         };
//         DynamicArr(const SELF_TYPE &src) : DynamicArr<ELEM_TYPE>(src.d, src.c) { this->copy(src); };
//         void swap(SELF_TYPE &rhs)
//         {
//             ArrVec<ELEM_TYPE>::swap(rhs);
//             using std::swap;
//             swap(this->n, rhs.n);
//             swap(this->c, rhs.c);
//         }
//         SELF_TYPE &operator=(const SELF_TYPE &rhs)
//         {
//             SELF_TYPE temp(rhs);
//             swap(temp);
//             return (*this);
//         }
//         ~DynamicArr() { this->release_array(this->n); };

//         void release_array()
//         {
//             if constexpr (!std::is_trivially_destructible<ELEM_TYPE>::value)
//                 for (auto i = 0; i < size; i++)
//                     this->v[i].~ELEM_TYPE();
//             std::free(this->v);
//             this->v = nullptr;
//         }

//         void release_entry_dynobj()
//         {
//             if constexpr (std::is_pointer<ELEM_TYPE>::value)
//                 for (auto i = 0; i < this->d; i++)
//                     if (this->v[i])
//                     {
//                         delete this->v[i];
//                         this->v[i] = nullptr;
//                     }
//         };

//         void release_entry_dynarr()
//         {
//             if constexpr (std::is_pointer<ELEM_TYPE>::value)
//                 for (auto i = 0; i < this->d; i++)
//                     if (this->v[i])
//                     {
//                         delete[] this->v[i];
//                         this->v[i] = nullptr;
//                     }
//         };

//         // DynamicArr DOES NOT support ELEM_TYPE being malloc pointers!!!

//         void Ones() { ArrVec<ELEM_TYPE>::Ones(n); };
//         void Zero() { ArrVec<ELEM_TYPE>::Zero(n); };

//         void push(const ELEM_TYPE &element)
//         {
//             // assert(!(n > d));
//             if (n == this->d)
//             {
//                 ELEM_TYPE *new_v = static_cast<ELEM_TYPE *>(std::malloc((this->vd + c) * sizeof(ELEM_TYPE)));
//                 if (!new_v)
//                     throw std::bad_alloc();

//                 if (this->v)
//                 {
//                     std::memcpy(new_v, this->v, this->d * sizeof(ELEM_TYPE));
//                     free(this->v); // IMPORTANT!!! NO deconstructor should be invoked here!!!
//                 }
//                 this->v = new_v;
//                 this->d = this->d + c;
//             }
//             this->v[n++] = element;
//         }

//         void fprintf(FILE *of, const char *prefix, INTE_TYPE length, const char *format) override
//         {
//             using std::fprintf;
//             INTE_TYPE dim = (length < n) ? length : n;

//             fprintf(of, "%s", prefix);
//             for (auto vec_ind = 0; vec_ind < dim; vec_ind++)
//                 fprintf(of, format, this->v[vec_ind]);
//             if (length < n)
//                 fprintf(of, "...\n\n");
//             else
//                 fprintf(of, "\n\n");
//         };

//         void printf(const char *prefix, INTE_TYPE length, const char *format) override
//         {
//             using std::printf;
//             INTE_TYPE dim = (length < n) ? length : n;
//             printf("%s", prefix);
//             for (auto vec_ind = 0; vec_ind < dim; vec_ind++)
//                 printf(format, this->v[vec_ind]);
//             if (length < n)
//                 printf("...\n\n");
//             else
//                 printf("\n\n");
//         };
//     };
// }
