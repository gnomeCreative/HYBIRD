#ifndef MACROS_H_
#define MACROS_H_

#include <cstdlib>

// a useful print info macro
#define INFO __FILE__ << ":" << __LINE__ << " "

// assert a condition
#define ASSERT(cond) if(!(cond)) { cout << INFO << "ERROR: Condition " #cond " not fulfilled!" << std::endl; exit(1); }
#define REQUEST(cond) if(!(cond)) { cout << INFO << "WARNING: Condition " #cond " not fulfilled!" << std::endl; }

// read a parameter from the config file, and if it is also given on the command line, overwrite it with that one
#define BASE_PARSE_CFG(cfg,type,var,name,def,def_typed) \
  if (!cfg.have_variable(name) && !commandLine.search("-" name))\
    cout << INFO << "WARNING: Parameter " name " is missing! Using default: " << def << std::endl;\
  type var = cfg(name, def_typed);\
  if (commandLine.search("-" name))\
    var = commandLine.next(var);

////// read a 3-vector parameter from the config file, and if it is also given on the command line, overwrite it with that one
/*
#define BASE_PARSE_CFG(cfg,type,var,name,def,def_typed) \
 if (!cfg.have_variable(name) && !commandLine.search("-" name))\
   cout << INFO << "WARNING: Parameter " name " is missing! Using default: " << def << std::endl;\
 type var = cfg(name, def_typed);\
 if (commandLine.search("-" name))\
   var = commandLine.next(var);
*/

//// read a member of vector parameter from the config file, and if it is also given on the command line, overwrite it with that one
#define BASE_PARSE_CFG_VEC(cfg,type,var,name,index,def,def_typed) \
  if (!cfg.have_variable(name))\
    cout << INFO << "WARNING: Parameter " name " is missing! Using default: " << def << std::endl;\
  if (commandLine.vector_variable_size("-" name)>0) { \
    var = commandLine("-" name, def_typed, index);} else { \
      type var = cfg(name, def_typed, index);\
  }
    


// some convenient parser abbreviations
#define MAKE_NAME_PARSE_CFG(cfg,type,var,name,def) BASE_PARSE_CFG(cfg,type,var,name,def,static_cast<type>(def))
#define MAKE_PARSE_CFG(cfg,type,var,def) MAKE_NAME_PARSE_CFG(cfg,type,var,#var,def)

// from Alessandro. In case the variable already exists
#define PARSE_CLASS_MEMBER(cfg,var,name,def) BASE_PARSE_CFG(cfg,{},var,name,def,def)
#define PARSE_CLASS_MEMBER_VEC(cfg,var,name,index,def) BASE_PARSE_CFG_VEC(cfg,{},var,name,index,def,def)

// parse an unused variable to prevent it from being added to the UFOs
#define DUMMY_PARSE_CFG(cfg,var) \
  cfg.have_variable(#var);\
  commandLine.search("-" #var);

// parse a 3-vector
#define MAKE_VECTOR_PARSE_CFG(cfg,var,def) \
  MAKE_NAME_PARSE_CFG(cfg, std::string, var##_expr, #var, #def "," #def "," #def);\
  mu::Parser p_##var;\
  p_##var.SetExpr(var##_expr);\
  int n_##var;\
  Real * var = p_##var.Eval(n_##var);\
  ASSERT(n_##var == 3)

#define MAKE_VECTOR_PARSE(var,def) MAKE_VECTOR_PARSE_CFG(shellCfgFile,var,def)

// parse a string variable and turn it to uppercase
#define STR_PARSE(var,def) \
  MAKE_NAME_PARSE_CFG(shellCfgFile, std::string, var##_name, #var, #def); \
  std::transform(var##_name.begin(), var##_name.end(), var##_name.begin(), ::toupper);

// a shortcut for assigning enums from strings
#define STR_TO_ENUM(var,val) if (var##_name.compare(#val) == 0) var = val;

// measure performance based on compile option (similar to libmesh_logging.h)
#ifdef MEASURE_PERFORMANCE
#include "perf_log.h"
#define PERFLOG_START(name) { perf_log.push(name); }
#define PERFLOG_STOP(name)  { perf_log.pop(name); }
#else
#define PERFLOG_START(name) {}
#define PERFLOG_STOP(name)  {}
#endif

// debugging helpers
#define HERE { cout << INFO << "Here!" << std::endl; }
#define EXIT { cout << INFO << "Exiting!" << std::endl; exit(0); }
#define PRINT(var) { cout << INFO << #var " = " << var << std::endl; }
#define CHECK(var) { if (!std::isfinite(var)) { PRINT(var) exit(1); } }

#endif /* MACROS_H_ */
