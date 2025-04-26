#ifndef JTF_MESH_CONF_H_
#define JTF_MESH_CONF_H_

#ifdef WIN32
#ifdef _EXPORTING
#define JTF_MESH_API __declspec(dllexport)
#else
#define JTF_MESH_API __declspec(dllimport)
#endif
#else
#  define JTF_MESH_API
#endif

#endif
