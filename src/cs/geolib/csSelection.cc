/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include <string>
#include <cstring>
#include <cstdio>
#include "csSelection.h"
#include "csSelectionField.h"
#include "csSelectionFieldDouble.h"
#include "csSelectionFieldInt.h"
#include "csException.h"
#include "geolib_string_utils.h"
#include "csFlexNumber.h"
#include "csVector.h"
#include "csGeolibUtils.h"
#include <cstdio>

using namespace cseis_geolib;

namespace {
  static const char HEADER_SEPARATOR = COLON;
  static const char SELECTION_SEPARATOR = FORWARD_SLASH;
  static const char FIELD_SEPARATOR = COMMA;
}

csSelection::csSelection( int numHeaders, type_t const* hdrTypes ) {
  myNumHeaders = 0;
  myHdrTypes   = NULL;
  mySelectionList = new csVector<csSelectionField const*>();
  myNumFieldsList = new csVector<int>();
  resetHeaders( numHeaders, hdrTypes );
}
csSelection::~csSelection() {
  if( mySelectionList != NULL ) {
    for( int i = 0; i < mySelectionList->size(); i++ ) {
      delete mySelectionList->at(i);
    }
    delete mySelectionList;
    mySelectionList = NULL;
  }
  if( myNumFieldsList != NULL ) {
    delete myNumFieldsList;
    myNumFieldsList = NULL;
  }
  if( myHdrTypes != NULL ) {
    delete [] myHdrTypes;
    myHdrTypes = NULL;
  }
}
//--------------------------------------------------------------------------
//
void csSelection::tokenize( std::string const& textIn, char separator, csVector<std::string>& tokenList ) {
  tokenList.clear();
  int counter = 0;
  int length = textIn.length();
  int pos1 = 0;

  while( counter < length ) {
    if( textIn[counter] == separator ) {
      tokenList.insertEnd( textIn.substr(pos1,counter-pos1) );
      counter++;
      pos1 = counter;
    }
    else {
      counter++;
    }
  }
  if( counter > pos1 ) {
    tokenList.insertEnd( textIn.substr(pos1,counter-pos1) );
  }
}
void csSelection::resetHeaders( int numHeaders, type_t const* hdrTypes ) {
  if( myHdrTypes != NULL ) {
    delete [] myHdrTypes;
    myHdrTypes = NULL;
  }
  myNumHeaders = numHeaders;
  myHdrTypes   = new type_t[myNumHeaders];
  memcpy( myHdrTypes, hdrTypes, myNumHeaders*sizeof(type_t) );

  //fprintf(stdout,"csSelection::resetHeaders: Header types: %d %s\n", myHdrTypes[0], cseis_geolib::csGeolibUtils::typeText(myHdrTypes[0]));
}
void csSelection::clear() {
  mySelectionList->clear();
  myNumFieldsList->clear();
}
//--------------------------------------------------------------------------------
void csSelection::add( std::string const& textIn ) {
  std::string text = replaceStr( textIn, " ", "" );

  csVector<std::string> fieldTokenList;
  csVector<std::string> headerTokenList;
  csVector<std::string> selectionTokenList;

  tokenize( text, SELECTION_SEPARATOR, selectionTokenList );
  int numSelections = selectionTokenList.size();
  if( numSelections < 1 ) throw csException("Syntax error in selection expression '%s': No valid selection found.", textIn.c_str());

  for( int iSelection = 0; iSelection < numSelections; iSelection++ ) {
    headerTokenList.clear();
    tokenize( selectionTokenList.at(iSelection), HEADER_SEPARATOR, headerTokenList );
    if( headerTokenList.size() != myNumHeaders ) {
      throw csException( "Syntax error in selection expression. Mismatch between number of headers and selections: %s", textIn.c_str() );
    }
    for( int iHeader = 0; iHeader < myNumHeaders; iHeader++ ) {
      if( headerTokenList.at(iHeader).length() == 0 ) {
        throw( csException("Syntax error in selection expression '%s': Header selection field #%d is empty (fields separated by character '%c').\nUse %c to select all header values", textIn.c_str(), iHeader+1, HEADER_SEPARATOR, STAR_CHAR) );
      }
      tokenize( headerTokenList.at(iHeader), FIELD_SEPARATOR, fieldTokenList );
      int numSelectionFields = fieldTokenList.size();
      myNumFieldsList->insertEnd(numSelectionFields);
      csSelectionField* selectionField;
      for( int iField = 0; iField < numSelectionFields; iField++ ) {
        if( myHdrTypes[iHeader] == TYPE_DOUBLE || myHdrTypes[iHeader] == TYPE_FLOAT ) {  // BUGFIX 29 Aug 2009: Treat FLOAT as DOUBLE (was: INT)
          selectionField = new csSelectionFieldDouble();
        }
        else {
          selectionField = new csSelectionFieldInt();
        }
        parse( fieldTokenList.at(iField), selectionField );
        mySelectionList->insertEnd( selectionField );
      }
    }
  }
}

//--------------------------------------------------------------------------------
bool csSelection::contains( csFlexNumber const* const values ) {
  int numSelectionFields = mySelectionList->size();
  int counterNumFields = 0;
  int counterSelectionFields = 0;
  bool isSelected = false;

  while( counterSelectionFields < numSelectionFields ) {
    bool isSelectedHeaders = true;
    for( int iHeader = 0; iHeader < myNumHeaders; iHeader++ ) {
      int numFields = myNumFieldsList->at(counterNumFields);
      counterNumFields += 1;
      bool isSelectedFields = false;
      for( int iField = 0; iField < numFields; iField++ ) {
        isSelectedFields = isSelectedFields || mySelectionList->at(counterSelectionFields)->contains( values[iHeader] );
        counterSelectionFields += 1;
      }
      isSelectedHeaders = isSelectedHeaders && isSelectedFields;  // All headers need to be selected
/*      if( !isSelectedFields ) {   // Every header must be selected true. Try next selection
        for( int ih = iHeader; ih < numHeaders; ih++ ) {  // Loop to set counters correctly
          int numFields = myNumFieldsList->at(counterNumFields++);
          counterSelectionFields += numFields;
        }
        break;
      } */
    }
    isSelected = isSelected || isSelectedHeaders;
  } // END while
    
  return isSelected;
}


//--------------------------------------------------------------------------------
void csSelection::dump() {
  int numSelectionFields = mySelectionList->size();
  int counterNumFields = 0;
  int counterSelectionFields = 0;

  int numSelections = mySelectionList[0].size();
  int counterSelections = 0;
  printf("==============================  (%d)\n", numSelections);
  while( counterSelectionFields < numSelectionFields ) {
    printf(" Selection: %d\n", counterSelections+1 );
    counterSelections += 1;
    for( int iHeader = 0; iHeader < myNumHeaders; iHeader++ ) {
      printf("    Header: %d\n", iHeader+1 );
      int numFields = myNumFieldsList->at(counterNumFields++);
      for( int iField = 0; iField < numFields; iField++ ) {
        printf("     Field: %d\n", iField+1 );
        mySelectionList->at(counterSelectionFields++)->dump();
      }
    }
  }

  printf(" END DUMP   ==============================\n");
}

//--------------------------------------------------------------------------------
//
void csSelection::parse( std::string const& textIn, csSelectionField* field ) {
  std::string text = replaceStr( textIn, " ", "" );

  csFlexNumber    value;
  csFlexNumber    rangeMin;
  csFlexNumber    rangeMax;
  csFlexNumber    rangeInc;
  csFlexNumber    selectionWidth;
  int op = csSelection::NONE;
  int selectType = csSelection::NONE;
  bool doInvert = false;

  csVector<csFlexNumber> valueList(3);
  int counter = 0;
  int pos1 = 0;
  int length = text.length();
  csFlexNumber valueTmp;

  while( counter < length ) {
    switch( text[counter] ) {
    case EXCL_CHAR:
      if( op != csSelection::NONE || selectType != csSelection::NONE || pos1 != counter )
        throw csSelectionException( std::string("Operator '!' in unexpected place"), textIn);
      if( doInvert ) throw csSelectionException("More than one exclamation mark ('NOT' operator !) found in selection.", textIn);
      doInvert = true;
      pos1 = counter+1;
      break;
    case GT_CHAR:
      if( op != csSelection::NONE || selectType != csSelection::NONE || pos1 != counter )
        throw csSelectionException("Operator '>' in unexpected place", textIn);
      if( counter < length-1 && text[counter+1] == EQ_CHAR ) {
        op = csSelection::OPERATOR_GREATER_EQUAL;
        counter++;
      }
      else {
        op = csSelection::OPERATOR_GREATER;
      }
      selectType = csSelection::SELECTION_OPERATOR;
      pos1 = counter+1;
      break;
    case LT_CHAR:
      if( op != csSelection::NONE || selectType != csSelection::NONE || pos1 != counter )
        throw csSelectionException("Operator '<' in unexpected place", textIn);
      if( counter < length-1 && text[counter+1] == EQ_CHAR ) {
        op = csSelection::OPERATOR_SMALLER_EQUAL;
        counter++;
      }
      else {
        op = csSelection::OPERATOR_SMALLER;
      }
      selectType = csSelection::SELECTION_OPERATOR;
      pos1 = counter+1;
      break;
    case DASH_CHAR:
      // Has to come directly behind a number to be the dash operator. Otherwise it has to be a number sign
      if( pos1 == counter ) break;  // This is surely a minus sign
      if( text[counter-1] == 'e' ) break;  // This dash is the sign of the 'e' notation in a number --> continue
      if( selectType != csSelection::NONE && op != csSelection::NONE )
        throw csSelectionException("Operator '-' in unexpected place", textIn);

      if( !valueTmp.convertToNumber( text.substr(pos1,counter-pos1) ) ) {
        throw csSelectionException( text.substr(pos1,counter-pos1), textIn );
      }
      valueList.insertEnd( valueTmp );

      selectType = csSelection::SELECTION_RANGE;
      pos1 = counter+1;
      break;
    case BRACKET_OPEN_CHAR:
      if( (selectType != csSelection::SELECTION_RANGE && op != csSelection::NONE) || pos1 == counter )
        throw csSelectionException("Bracket '(' in unexpected place", textIn);
      selectType = csSelection::SELECTION_RANGE_INC;

      if( !valueTmp.convertToNumber( text.substr(pos1,counter-pos1) ) ) {
        throw csFlexNumberException( text.substr(pos1,counter-pos1) );
      }
      valueList.insertEnd( valueTmp );

      pos1 = counter+1;
      break;
    case BRACKET_CLOSE_CHAR:
      if( (selectType != csSelection::SELECTION_RANGE_INC && op != csSelection::NONE) || pos1 == counter )
        throw csSelectionException("Bracket ')' in unexpected place", textIn);

      selectType = csSelection::SELECTION_RANGE_INC;

      if( !valueTmp.convertToNumber( text.substr(pos1,counter-pos1) ) ) {
        throw csFlexNumberException( text.substr(pos1,counter-pos1) );
      }
      valueList.insertEnd( valueTmp );

      pos1 = counter+1;
      break;
    case PLUS_CHAR:
      if( (selectType != csSelection::SELECTION_RANGE_INC && selectType != csSelection::NONE) || op != csSelection::NONE )
        throw csSelectionException("Selection width operator '+' found in unexpected place.", textIn);
      if( selectType == csSelection::NONE ) {
        selectType = csSelection::SELECTION_SINGLE;
        if( valueTmp.convertToNumber( text.substr(pos1,counter-pos1) ) ) {
          valueList.insertEnd( valueTmp );
        }
        else {
          throw csFlexNumberException( text.substr(pos1,counter-pos1) );
        }
      }
      if( !valueTmp.convertToNumber( text.substr(counter,length-counter) ) ) {
        throw csFlexNumberException( text.substr(counter,length-counter) );
      }
      valueList.insertEnd( valueTmp );
      counter = length-1;
      pos1    = counter+1;
      break;
    case STAR_CHAR:
      selectType = csSelection::SELECTION_ALL;
      pos1 = counter+1;
      break;
    } // End switch
    counter++;
  } // End while

  if( pos1 != counter ) {
    if( valueTmp.convertToNumber( text.substr(pos1,counter-pos1) ) ) {
      valueList.insertEnd( valueTmp );
      if( selectType == csSelection::NONE ) {
        selectType = csSelection::SELECTION_SINGLE;
      }
    }
    else {
      throw csFlexNumberException( text.substr(pos1,counter-pos1) );
    }
  }

  switch( selectType ) {
  case csSelection::SELECTION_SINGLE:
    if( valueList.size() == 1 ) {
      field->setValue( valueList.at(0) );
    }
    else {
      selectType = csSelection::SELECTION_RANGE;
      field->setRange( valueList.at(0) - valueList.at(1), valueList.at(0) + valueList.at(1) );
    }
    break;
  case csSelection::SELECTION_OPERATOR:
    if( valueList.size() != 1 )
      throw csSelectionException("Single selection expects one operand.", textIn);
    field->setValue( valueList.at(0) );
    break;
  case csSelection::SELECTION_RANGE:
    if( valueList.size() != 2 )
      throw csSelectionException("Range selection expects two operands.", textIn);
    if( valueList.at(0) < valueList.at(1) ) {
      field->setRange( valueList.at(0), valueList.at(1) );
    }
    else {
      field->setRange( valueList.at(1), valueList.at(0) );
    }
    break;
  case csSelection::SELECTION_RANGE_INC:
    if( valueList.size() < 3 )
      throw csSelectionException("Range selection with increment expects three operands (plus one for width).", textIn);
    if( valueList.at(2) == 0.0 ) throw csSelectionException("Range selection increment cannot be zero.", textIn);
    if( (valueList.at(1) - valueList.at(0)) / valueList.at(2) < 0.0 ) {
      throw csSelectionException("Range selection increment and min/max range inconsistent.", textIn);
    }
    if( valueList.at(0) < valueList.at(1) ) {
      field->setRange( valueList.at(0), valueList.at(1), valueList.at(2) );
    }
    else {
      field->setRange( valueList.at(1), valueList.at(0), valueList.at(2)*-1.0 );
    }
    if( valueList.size() > 3 ) {
      selectType = csSelection::SELECTION_RANGE_INC_WIDTH;
      field->setWidth( valueList.at(3) );
    }
    break;
  case csSelection::SELECTION_ALL:
    break; // Nothing
  default:
    throw csSelectionException("Program bug.", textIn);
  }  // END switch

  field->set( selectType, op, doInvert );
//  printf("Select type is...%d\n", selectType);
}
//--------------------------------------------------------------------
//
csSelectionException::csSelectionException( std::string const& message, std::string const& expression ) : csException() {
  myMessage = std::string("Syntax error in selection expression '") + expression + std::string("':\n") + message;
}
csSelectionException::~csSelectionException() {
}

