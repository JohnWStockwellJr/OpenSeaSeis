/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_USER_CONSTANT_H
#define CS_USER_CONSTANT_H

#include <string>

namespace cseis_system {

/**
* User constant
* A user constant is a constant supplied by the user in the job flow, using a &define statement
* A user constant is inherently a string that then replaces all occurences of the
* user constant name, surrounded by &, found in the job flow
*
* @author Bjorn Olofsson
* @date   2007
*/
class csUserConstant {
public:
  static char const LETTER_DEFINE = '&';
public:
  csUserConstant();
  csUserConstant( std::string const& name, std::string const& value );
  csUserConstant( csUserConstant const& obj );
  ~csUserConstant();
  csUserConstant& operator=( csUserConstant const& obj );
  bool operator==( csUserConstant const& obj );
  std::string pureName();
  static std::string const defineWord() {
    static std::string const WORD_DEFINE = std::string("&define");
    return WORD_DEFINE;
  }
  static std::string const tableWord() {
    static std::string const TABLE_DEFINE = std::string("&table");
    return TABLE_DEFINE;
  }
  void setValue( std::string theValue );
  std::string name;
  int nameLength;
  std::string value;
};

} // namespace

#endif

