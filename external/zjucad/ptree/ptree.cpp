#include <iostream>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/info_parser.hpp>

#include "ptree.h"

namespace zjucad
{
	using namespace boost;

	static const char *error_info_brief = 
    "argument error: prog or tet haven't been set, please use A=B format\n"
    " prog ""                  "
    "<frame, frame2vtk, cut_tet, int_param, find_singularities,find_singularities2, int_pts>\n"
    " tet ""                  "
    "file name to the input tet model";

	int load(const std::string & file_path, ptree & pt)
	{
    size_t last_dot_pos = file_path.find_last_of(".");
	  
	  if(last_dot_pos == std::string::npos)
      throw std::logic_error("argument error: do not support this format, only support [xml/json/ini/info]");
	  
	  std::string file_format = file_path.substr(last_dot_pos + 1);
	  
	  if(file_format == "xml" || file_format == "XML")
      read_xml(file_path, pt);
	  else if(file_format == "json" || file_format == "JSON")
      read_json(file_path, pt);
	  else if(file_format == "ini" || file_format == "INI")
      read_ini(file_path, pt);
	  else if(file_format == "info" || file_format == "INFO")
      read_info(file_path, pt);
	  else 
      throw std::logic_error("argument error: do not support this format, only support [xml/json/ini/info]");
	  return 0;
	}

    int load(const std::wstring & file_path, wptree & pt)
    {
        throw std::logic_error("can not open config file");
        return __LINE__;

//    size_t last_dot_pos = file_path.find_last_of(L".");

//      if(last_dot_pos == std::wstring::npos)
//      throw std::logic_error("argument error: do not support this format, only support [xml/json/ini/info]");

//      std::wstring file_format = file_path.substr(last_dot_pos + 1);

//      std::wifstream ifs(file_path.c_str());
//      if(file_format == L"xml" || file_format == L"XML")
//      read_xml(ifs, pt);
//      else if(file_format == L"json" || file_format == L"JSON")
//      read_json(ifs, pt);
//      else if(file_format == L"ini" || file_format == L"INI")
//      read_ini(ifs, pt);
//      else if(file_format == L"info" || file_format == L"INFO")
//      read_info(ifs, pt);
//      else
//      throw std::logic_error("argument error: do not support this format, only support [xml/json/ini/info]");
//      return 0;
    }
	
	int save(const std::string & file_path, const ptree & pt)
	{
    size_t last_dot_pos = file_path.find_last_of(".");
	  
	  if(last_dot_pos == std::string::npos)
      throw std::logic_error("argument error: do not support this format, only support [xml/json/ini/info]");
	  
	  std::string file_format = file_path.substr(last_dot_pos + 1);
  
	  if(file_format == "xml" || file_format == "XML")
      write_xml(file_path, pt);
	  else if(file_format == "json" || file_format == "JSON")
      write_json(file_path, pt);
	  else if(file_format == "ini" || file_format == "INI")
      write_ini(file_path, pt);
	  else if(file_format == "info" || file_format == "INFO")
      write_info(file_path, pt);
	  else 
      throw std::logic_error("argument error: do not support this format, only support [xml/json/ini/info]");
	  return 0;
	}

    int save(const std::wstring & file_path, const wptree & pt)
    {
        throw std::logic_error("can not open config file");
        return __LINE__;
//    size_t last_dot_pos = file_path.find_last_of(L".");

//      if(last_dot_pos == std::wstring::npos)
//      throw std::logic_error("argument error: do not support this format, only support [xml/json/ini/info]");

//      std::wstring file_format = file_path.substr(last_dot_pos + 1);

//      std::wofstream ofs(file_path.c_str());
//      if(file_format == L"xml" || file_format == L"XML")
//      write_xml(ofs, pt);
//      else if(file_format == L"json" || file_format == L"JSON")
//      write_json(ofs, pt);
//      else if(file_format == L"ini" || file_format == L"INI")
//      write_ini(ofs, pt);
//      else if(file_format == L"info" || file_format == L"INFO")
//      write_info(ofs, pt);
//      else
//      throw std::logic_error("argument error: do not support this format, only support [xml/json/ini/info]");
//      return 0;
    }
	
	int read_cmdline(int argc, char *argv[], ptree &pt)
	{
	  const size_t token_len = strlen("=");
	  bool need_help = false;
	  /// first to find the "config=xxx" term
	  for(int i = 0; i < argc; ++i)
      {
        std::string k_v = argv[i];
        size_t pos = k_v.find("=");
        if(pos != std::string::npos && k_v.substr(0, pos) == "config")
          {
            std::string value = k_v.substr(pos + token_len);
            load(value,pt);
          }
        if(k_v == "help" || k_v == "--help" || k_v == "-h")
          {
            need_help = true;
            break;
          }
        /// To ensure this application can get xml format data std::cin, there is a trick.
        /// I make a fake input data to gather sereral input section into one section, 
        /// and replace the cin.rdbuf with this fake data.
        if(k_v == "cin" || k_v == "--cin" || k_v == "-c")
          {
            std::string temp;
            while(std::cin.good())
              {
                std::string str;
                getline(std::cin,str);
                temp += str;
              }
            std::cin.clear();
            std::istringstream oss(temp.c_str());
            std::cin.rdbuf(oss.rdbuf());
			  
            read_xml(std::cin,pt);
          }
      }
	  
	  for(int i = 0; i < argc; ++i)
      {
        std::string k_v = argv[i];
        size_t pos = k_v.find("=");
        if(pos != std::string::npos)
          {
            std::string path = k_v.substr(0, pos);
			  
            if(path == "config") 
              continue;
			  
            std::string value = k_v.substr(pos + token_len);
            // default path is key.value
            const size_t last_dot = path.find_last_of(".");
            if(last_dot == std::string::npos)
              {
                pt.put(path + ".value", value);
                pt.put(path + ".is_required","y");
              }else{
                std::string last_term = path.substr(last_dot+1, path.size());
                if(last_term != "value" && last_term != "is_required"){
                    pt.put(path + ".value",value);
                    pt.put(path + ".is_required", "y");
                }else pt.put(path, value);
            }
          }
      }
	  // to add help info into the ptree, so that, user can get help message when he meets an excption
	  //init_help_info(pt);
	  if(need_help)
      throw std::logic_error("Help:");
	  return 0;
	}
	
	void show_usage_info(std::ostream & out, const ptree &pt)
	{
	  const size_t total_space_num = 24;
	  const std::string cin_info = "-c [--cin]";
	  const std::string space_after_cin((total_space_num - cin_info.length()) > 0?total_space_num - cin_info.length():0,' ');
	  out << cin_info << space_after_cin << "accept the input xml data from std::cin" << std::endl;
	  
	  for (ptree::const_iterator it = pt.begin(); it != pt.end(); ++it)
      {
        std::string desc_path = it->first.data();
        std::ostringstream oss;
        if(has(desc_path + ".desc",pt))
          {
            if(   pt.get<std::string>(desc_path + ".is_required","y") == "n" 
                  || pt.get<std::string>(desc_path + ".is_required","y") == "N" )
              oss << '[' << it->first.data() << ']';
            else
              oss << ' ' << it->first.data() << ' ' ;
            out << oss.str();
            if(oss.tellp() < total_space_num)
              {
                std::string space(total_space_num - oss.tellp(), ' ');
                out << space ;
              }
            out << pt.get<std::string>(desc_path + ".desc","") << std::endl;
          }
      }
	}

    void show_usage_info(std::wostream & out, const wptree &pt)
    {
      const size_t total_space_num = 24;
      const std::wstring cin_info = L"-c [--cin]";
      const std::wstring space_after_cin((total_space_num - cin_info.length()) > 0?total_space_num - cin_info.length():0,' ');
      out << cin_info << space_after_cin << L"accept the input xml data from std::cin" << std::endl;

      for (wptree::const_iterator it = pt.begin(); it != pt.end(); ++it)
      {
        std::wstring desc_path = it->first.data();
        std::wostringstream oss;
        if(has(desc_path + L".desc",pt))
          {
            if(   pt.get<std::wstring>(desc_path + L".is_required",L"y") == L"n"
                  || pt.get<std::wstring>(desc_path + L".is_required",L"y") == L"N" )
              oss << L'[' << it->first.data() << L']';
            else
              oss << L' ' << it->first.data() << L' ' ;
            out << oss.str();
            if(oss.tellp() < total_space_num)
              {
                std::wstring space(total_space_num - oss.tellp(), L' ');
                out << space ;
              }
            out << pt.get<std::wstring>(desc_path + L".desc",L"") << std::endl;
          }
      }
    }

	bool has(const std::string& name, const ptree &pt)
	{
    // is pt.get_child_optional("name") better?
		std::string value = pt.get<std::string>(name,"");
		if(value != "")
			return true;
		else return false;
	}

    bool has(const std::wstring& name, const wptree &pt)
    {
    // is pt.get_child_optional("name") better?
        std::wstring value = pt.get<std::wstring>(name,L"");
        if(value != L"")
            return true;
        else return false;
    }
}
