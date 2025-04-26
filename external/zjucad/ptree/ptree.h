#ifndef _ZJUCAD_PTREE_H_
#define _ZJUCAD_PTREE_H_

/**
 * Author: Tengfei Jiang, based on hjlib::arg_opt
 * Notice:
 *        1. The path in property_tree is separated by dot "."         
 *        2. In my desgin, each entry should have three subentries, "value", "desc", "is_required". 
 *                 "value" stands for the key-value; 
 *                  "desc" stands for the description; 
 *           "is_required" stands for the key-priority, you should only choose from [y/n], default is y. 
 *        3. Only support one level list: (need to be extended) 
 *         --------------------------------------------------------> 
 *            A |            B |            C |              D | 
 *         val des  is    val des is     val des is       val des is
 */

#include <boost/property_tree/ptree.hpp>
#include <string>

#ifdef PTREE_COMMON_EXPORT
#define PTREE_COMMON_API __declspec(dllexport)
#else
#define PTREE_COMMON_API
#endif

namespace zjucad
{
  using boost::property_tree::ptree;
  using boost::property_tree::wptree;
  /** 
   * Read configure file into property_tree, only support
   * "xml/json/ini/info" format. Do not suppot full features of above
   * formats. Check boost/property_tree for details.
   *
   * @param file_path the pah of file 
   * @param pt ptree
   * 
   * @return sucess return 0, others throw an exception
   */
  PTREE_COMMON_API
  int load(const std::string & file_path, ptree &pt);

  PTREE_COMMON_API
  int load(const std::wstring & file_path, wptree &pt);

  /** 
   * Write out configure file from property_tree, only support
   * "xml/json/ini/info" format. Do not suppot full features of above
   * formats. Check boost/property_tree for details.
   *
   * @param file_path the pah of file 
   * @param pt ptree
   * 
   * @return sucess return 0, others throw an exception
   */
  PTREE_COMMON_API 
  int save(const std::string & file_path, const ptree &pt);

  PTREE_COMMON_API
  int save(const std::wstring & file_path, const wptree &pt);

  /** 
   * Read parameters from cmdline, only support "A=B" format. "A=B" or
   * "A.value=B" will add entry <A.value,B> into property tree
   * "A.desc=B" will add entry <A.desc,B> into property tree If find
   * "config=xxx", that means there is a config file which need to be
   * loaded first, but with the lower priority.  If there are other
   * inputs from cmdline at the same time, the input has higher
   * priority. That means the cmdline input will correct/add entries
   * within the configure file.
   *
   * @param argc number of main parameters
   * @param argv main parameters
   * @param pt ptree
   * 
   * @return sucess return 0, others throw an exception
   */
  PTREE_COMMON_API 
  int read_cmdline(int argc, char *argv[], ptree &pt);

			
  /** 
   * Add description and is_required info into the ptree, to show help
   * info when an excption was caught Notice that, the build-in info
   * has the lowest priority, if the info has been input through
   * cmdline or configure file, the entry will not be modified.  Only
   * the empty or new item will be added.  @param pt ptree
   */
  PTREE_COMMON_API 
  void show_usage_info(std::ostream & out, const ptree &pt);

  PTREE_COMMON_API
  void show_usage_info(std::wostream & out, const wptree &pt);
				
  PTREE_COMMON_API
  bool has(const std::string& name, const ptree &pt);

  PTREE_COMMON_API
  bool has(const std::wstring& name, const wptree &pt);
}

#endif
