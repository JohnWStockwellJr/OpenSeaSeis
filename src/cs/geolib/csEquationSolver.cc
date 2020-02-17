/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#include "csEquationSolver.h"
#include "csMathConstant.h"
#include "geolib_math.h"
#include "geolib_string_utils.h"
#include "csMathFunction.h"
#include "csStack.h"
#include "csVector.h"

using namespace cseis_geolib;

//------------------------------------------------------------
//
csEquationSolver::csEquationSolver( std::string nameSpecialFloatArray, int numValuesSpecialFloatArray ) {
  init();
  myNumValuesSpecialFloatArray = numValuesSpecialFloatArray;
  myNameSpecialFloatArray = nameSpecialFloatArray;
  myIsSpecialFloatArray = true;
}
csEquationSolver::csEquationSolver() {
  init();
}
void csEquationSolver::init() {
  myIsSpecialFloatArray   = false;
  myNameSpecialFloatArray = "";
  myNumValuesSpecialFloatArray = 0;
  myNumUserConstants = 0;
  myUserConstants = NULL;
  myUserConstantNames = NULL;
  myValueStack = NULL;
  myExpressionList = new csVector< csVector<csToken> >;
}
//------------------------------------------------------------
//
csEquationSolver::~csEquationSolver() {
  myNumUserConstants = 0;
  if( myUserConstants ) {
    delete [] myUserConstants;
    myUserConstants = NULL;
  }
  if( myUserConstantNames ) {
    delete [] myUserConstantNames;
    myUserConstantNames = NULL;
  }
  if( myValueStack ) {
    delete [] myValueStack;
    myValueStack = NULL;
  }
  if( myExpressionList ) {
    delete myExpressionList;
    myExpressionList = NULL;
  }
}

//------------------------------------------------------------
//
std::string const& csEquationSolver::getErrorMessage() const {
  return myErrorMessage;
}

//------------------------------------------------------------
//
bool csEquationSolver::prepare( std::string const& equation ) {
  return prepare( equation, NULL, 0, NULL );
}
bool csEquationSolver::prepare( std::string const& equation, std::string const* userConstantNames, int nConstants ) {
  return prepare( equation, userConstantNames, nConstants, NULL );
}
bool csEquationSolver::prepare( std::string const& equation, std::string const* userConstantNames, int nConstants, double const* userConstants ) {
  std::string expression = removeChar( toLowerCase( equation ), ' ' );

  // Preliminary setup of user constant names. These lists will be reduced to the constants that are actually referenced in the equation
  myNumUserConstants  = nConstants;
  myUserConstants     = new double[myNumUserConstants];
  myUserConstantNames = new std::string[myNumUserConstants];
  for( int i = 0; i < myNumUserConstants; i++ ) {
    myUserConstants[i] = 0.0;
    myUserConstantNames[i] = toLowerCase( userConstantNames[i] );
  }

  if( userConstants != NULL ) {
    for( int i = 0; i < myNumUserConstants; i++ ) {
      myUserConstants[i] = userConstants[i];
    }
  }

  try {
    csVector<csToken> tokenList;
    tokenizeExpression( expression, tokenList );
    prepareExpressionList( tokenList );
  }
  catch( ExpressionException& e ) {
    //    string message = std::string("Unable to evaluate equation '") + equation + std::string("':\n ");
    //   myErrorMessage = message + string(e.getMessage());
    myErrorMessage = e.getMessage();
    return false;
  }
  catch( ConstantException& e ) {
    if( DEBUG_EQ_SOLVER ) printf("ConstantException caught... \n");
    myErrorMessage = std::string("Program bug: ") + e.getMessage() + equation;
    return false;
  }
  catch( EquationException& e ) {
    //if( DEBUG_EQ_SOLVER )
    myErrorMessage = e.getMessage() +  std::string("\nEquation: ") + equation;
    return false;
  }
  catch(...) {
    if( DEBUG_EQ_SOLVER ) printf("Unknown Exception caught... \n");
    myErrorMessage = std::string("Program bug: Unknown exception caught.\nEquation: ") + equation;
    return false;
  }
  return true;
}
//------------------------------------------------------------
//
csMathFunction const* csEquationSolver::getMathFunction( std::string functionName ) {
  csMathFunction const* const presets = csMathFunction::presets();
  for( int i = 0; i < csMathFunction::NUM_PRESETS; i++ ) {
    if( !functionName.compare( presets[i].name ) ) {
      return &presets[i];
    }
  }
  return NULL;
}
csMathConstant const* csEquationSolver::getMathConstant( std::string constantName ) {
  csMathConstant const* const presets = csMathConstant::presets();
  for( int i = 0; i < csMathConstant::NUM_PRESETS; i++ ) {
    if( !constantName.compare( presets[i].name ) ) {
      return &presets[i];
    }
  }
  return NULL;
}

//------------------------------------------------------------
//
int csEquationSolver::getUserConstantIndex( std::string const& identifierName ) const {
  if( myUserConstantNames != NULL ) {
    for( int i = 0; i < myNumUserConstants; i++ ) {
      if( !myUserConstantNames[i].compare( identifierName ) ) {
        return i;
      }
    }
  }
  printf("Constant exception thrown on '%s'...\n", identifierName.c_str());
  throw( ConstantException( std::string("Constant not found: '") + identifierName) );
}

//------------------------------------------------------------
//
bool csEquationSolver::isUserConstant( std::string const& identifierName ) const {
  if( myUserConstantNames != NULL ) {
    for( int i = 0; i < myNumUserConstants; i++ ) {
      if( !myUserConstantNames[i].compare( identifierName ) ) {
        return true;
      }
    }
  }
  return false;
}

//------------------------------------------------------------
//
bool csEquationSolver::isOperator(char c) const {
  if( c == '+' || c == '-' || c == '*' || c == '/' )
    return true;
  else
    return false;
}

//------------------------------------------------------------
//
bool csEquationSolver::isPrefixOperator( std::string const& expression, int operatorPos ) const {
  if( operatorPos == 0 ) {
    return true;
  }
  else if( operatorPos == (int)expression.length()-1 ) {
    return false;
  }
  else {
    char cPrev = expression[operatorPos-1];
    if( isOperator(cPrev) || cPrev == '(' || cPrev == ',' ) {
      return true;
    }
    else {
      return false;
    }
  }
}

//------------------------------------------------------------
//
int csEquationSolver::getCharType( char c ) const {
  if( isDigit(c) ) {
    return DIGIT;
  }
  else if( c == '(' || c == ')' ) {
    return BRACKET;
  }
  else if( isOperator(c) ) {
    return OPERATOR;
  }
  else if( isLetter(c) || c == '_' ) {
    return LETTER;
  }
  else if( c == ',' ) {
    return COMMA;
  }
  else if( c == '.' ) {
    return DOT;
  }
  else {
    return UNKNOWN;
  }
}

//------------------------------------------------------------
//
void csEquationSolver::tokenizeExpression( std::string expression, csVector<csToken>& tokenList )  {
  int fullLength = expression.length();
  int counter = 0;
  int numCloseBrackets = 0;
  int numOpenBrackets = 0;
  char c;
  
  while( counter < fullLength ) {
    if( DEBUG_EQ_SOLVER ) printf("Processing %d of %d\n", counter, fullLength );
    c = expression[ counter ];
    int charType = getCharType( c );
    
    switch( charType ) {
    case BRACKET:
      if( DEBUG_EQ_SOLVER ) printf( "Bracket type: '%c'\n", c );
      if( c == '(' ) {
        tokenList.insertEnd( csToken( c, BRACKET_OPEN, counter ) );
        numOpenBrackets++;
      }
      else {
        tokenList.insertEnd( csToken( c, BRACKET_CLOSE, counter ) );
        numCloseBrackets++;
      }
      counter += 1;
      break;
      //--------------------------
    case OPERATOR:
      if( DEBUG_EQ_SOLVER ) printf( "Operator type: '%c'\n", c );
      if( c == '+' || c == '-' ) {
        if( isPrefixOperator( expression, counter ) ) {
          // Replace prefix - with -1.0 *
          if( c == '-' ) {
            tokenList.insertEnd( csToken( '-', OPERATOR_PREFIX_MIN, counter, counter ) );
          }
          else {
            throw(  ExpressionException("Syntax error. Prefix '+' operator not supported",
                                        csToken( " ", UNKNOWN, counter, counter+1 ) ) );
          }
        }
        else {
          tokenList.insertEnd( csToken( c, OPERATOR_PLUS_MINUS, counter ) );
        }
      }
      else {
        tokenList.insertEnd( csToken( c, OPERATOR_MULT_DIV, counter ) );
      }
      counter += 1;
      break;
      //--------------------------
    case DOT:
    case DIGIT:
      {
        if( DEBUG_EQ_SOLVER ) printf( "Digit/Dot type: '%c'\n", c );
        int beginIndex = counter;

        // 1. Number starts with a series of digits:
        while( counter < fullLength && isDigit( expression[ counter ] ) ) {
          counter++;
        }

        // 2. optional dot for floating point number
        if( counter < fullLength && expression[ counter ] == '.' ) {
          counter++;
          while( counter < fullLength && isDigit( expression[ counter ] ) ) {
            counter++;
          }
        }
      
        // 3. Optional 'e' notation
        if( counter < fullLength && expression[ counter ] == 'e' ) {
          counter++;
          if( counter >= fullLength || ( expression[ counter ] != '+' && expression[ counter ] != '-' ) ) {
            throw(  ExpressionException( "Wrong number format (missing sign in 'e' notation)",
                                         csToken( 0.0, NUMBER, beginIndex, counter) ) );
          }
          counter++;
          int checkIndex = counter;
          while( counter < fullLength && isDigit( expression[ counter ] ) ) {
            counter++;
          }
          if( counter == checkIndex ) {
            throw( ExpressionException( "Wrong number format (no number in 'e' notation)",
                                        csToken(0.0, NUMBER, beginIndex, counter) ) );
          }
        }
        if( DEBUG_EQ_SOLVER ) printf( "    Number index: %d %d\n", beginIndex, counter );
        double newNumber = atof( expression.substr( beginIndex, counter ).c_str() );
        tokenList.insertEnd( csToken( newNumber, NUMBER, beginIndex, counter ) );
      }
      break;
      //--------------------------
    case LETTER:
      {
        if( DEBUG_EQ_SOLVER ) printf( "Letter type: '%c'\n", c );
        int subCounter = counter + 1;
        while( subCounter < fullLength && (isLetter( expression[ subCounter ] ) || expression[ subCounter ] == '.' ||
                                           isDigit(expression[ subCounter ])) ) {
          subCounter++;
        }
        std::string newWord = expression.substr( counter, subCounter-counter );
        csMathFunction const* mathFunction;
        csMathConstant const* mathConstant;
        if( DEBUG_EQ_SOLVER ) printf( "   Word is: '%s', counters: %d %d\n", newWord.c_str(), counter, subCounter );
        // ...is a variable:
        if( isUserConstant( newWord ) ) {
          if( DEBUG_EQ_SOLVER ) printf( "new word is user constant! %s\n", newWord.c_str() );
          tokenList.insertEnd( csToken( newWord, USER_CONSTANT, getUserConstantIndex(newWord), counter, subCounter ) );
        }
        // ...is a math function:
        else if( (mathFunction = getMathFunction(newWord)) != NULL ) {
          if( DEBUG_EQ_SOLVER ) printf( "   csMath function: %s %d\n", mathFunction->name.c_str(), mathFunction->key );
          tokenList.insertEnd( csToken( mathFunction, counter, subCounter ) );
        }
        else if( (mathConstant = getMathConstant(newWord)) != NULL ) {
          if( DEBUG_EQ_SOLVER ) printf( "   csMath constant: %s %d\n", mathConstant->name.c_str(), mathConstant->key );
          tokenList.insertEnd( csToken( mathConstant->value, NUMBER, counter, subCounter ) );
        }
        else {
          if( DEBUG_EQ_SOLVER ) printf( "Unknown identifier: %s\n", newWord.c_str() );
          throw( ExpressionException( "Unknown identifier: '" + newWord + "'",
                                      csToken( " ", UNKNOWN, counter, subCounter ) ) );
        }
        counter = subCounter;
      }
      break;
      //--------------------------
    case COMMA:
      //if( DEBUG_EQ_SOLVER ) System.out.println( "Comma type: " + c );
      //    throw( new ExpressionException("Commas/multiple argument functions are currently not supported yet",
      //        new csTokenObject(',',UNKNOWN,counter)) );
      tokenList.insertEnd( csToken( ',', COMMA_OPERATOR, counter ) );
      counter += 1;
      break;
      //--------------------------
    case UNKNOWN:
      if( DEBUG_EQ_SOLVER ) printf("Unknown character: '%c'\n", c );
      char c2[4];
      c2[0] = '"';
      c2[1] = c;
      c2[2] = '"';
      c2[3] = '\0';
      throw( ExpressionException( std::string("Unknown letter in equation: ") + std::string(c2), csToken(c, UNKNOWN, counter)) );
    }
  }

  if( numCloseBrackets != numOpenBrackets ) {
    throw( ExpressionException( "Bracket mismatch. Open brackets/closed brackets: " ) ); // + numOpenBrackets +
    //     ", closed brackets: " + numCloseBrackets ) );
  }
  if( DEBUG_EQ_SOLVER ) {
    printf( "---------------------------------\n" );
    printf( " csTokenized queue/stack: \n" );
    for( int i = 0; i < tokenList.size(); i++ ) {
      csToken t = tokenList.at(i);
      if( t.typeValue == DOUBLE ) {
        printf( "Double: %f\n", t.valDouble );
      }
      else if( t.typeValue == STRING ) {
        printf( "String: %s\n", t.valString.c_str() );
      }
      else if( t.typeValue == CHAR ) {
        printf( "Char:   '%c'\n", t.valChar );
      }
      else if( t.typeValue == FUNC ) {
        printf( "Func:   '%s'  %d\n", t.function->name.c_str(), t.function->key );
      }
      else {
        printf( "????:   \n" );
      }
    }
  }
}
//--------------------------------------------------------------------------------
//
void csEquationSolver::prepareUserConstants( csVector<std::string>* const constantList ) {
  int nExpressions = myExpressionList->size();
  csToken* tokenPtr;
  for( int iex = 0; iex < nExpressions; iex++ ) {
    int nTokens = myExpressionList->at(iex).size();
    for( int itoken = 0; itoken < nTokens; itoken++ ) {
      tokenPtr = &myExpressionList->at(iex).at(itoken);
      if( tokenPtr->type == USER_CONSTANT ) {
        tokenPtr->valInt = constantList->size();  // TEMP!!! !CHANGE! just for testing...
        constantList->insertEnd( tokenPtr->valString );
      }
    }
  }

  if( myUserConstants ) delete [] myUserConstants;
  if( myUserConstantNames ) delete [] myUserConstantNames;

  myNumUserConstants = constantList->size();
  myUserConstants = new double[myNumUserConstants];
  myUserConstantNames = new std::string[myNumUserConstants];

  for( int i = 0; i < myNumUserConstants; i++ ) {
    myUserConstants[i] = 0.0;
    myUserConstantNames[i] = toLowerCase( constantList->at(i) );
  }
}

//--------------------------------------------------------------------------------
//
void csEquationSolver::setUserConstants( double const* userConstants, int nConstants ) {
  if( nConstants != myNumUserConstants ) {
    //    fprintf(stderr,"No of constants: %d\n\n", myNumUserConstants);
  }
  else {
    for( int i = 0; i < myNumUserConstants; i++ ) {
      myUserConstants[i] = userConstants[i];
    }
  }
}

//--------------------------------------------------------------------------------
//
double csEquationSolver::solve() {
  int size = myExpressionList->size();
  //  double result;
  csToken* tokenArg1Ptr;
  csToken* tokenOp1Ptr;
  csToken* tokenOp2Ptr;
  double arg1, arg2, arg3;
  csToken token;
  csToken tokenArgument;
  tokenArgument.type = NUMBER;

  for( int iex = 0; iex < size; iex++ ) {
    int nTokens = myExpressionList->at(iex).size();
    tokenArg1Ptr = &myExpressionList->at(iex).at(0);
    //    printf("Solve expression: #%d, nTokens: %d\n", iex, nTokens);

    if( tokenArg1Ptr->type == FUNCTION ) {
      arg1 = getValue( &(myExpressionList->at(iex).at(1)) );
      if( nTokens == 2 ) {
        //  printf("  input function '%s' %f\n", tokenArg1Ptr->function->name.c_str(), arg1 );
        arg1 = tokenArg1Ptr->function->method1ArgPtr( arg1 );
      }
      
      else {
        arg2 = getValue( &(myExpressionList->at(iex).at(2)) );
        //  printf("  input function '%s' %f %f\n", tokenArg1Ptr->function->name.c_str(), arg1, arg2 );
        arg1 = tokenArg1Ptr->function->method2ArgPtr( arg1, arg2 );
      }
      //      printf("  function '%s' result: %f\n", tokenArg1Ptr->function->name.c_str(), arg1 );
    }
    else {
      arg1 = getValue( tokenArg1Ptr );
      //      printf("  arg1: %f\n", arg1);

      if( nTokens > 1 ) {
        int itoken = 1;
        tokenOp1Ptr = &myExpressionList->at(iex).at(itoken++);
        arg2 = getValue( &(myExpressionList->at(iex).at(itoken++)) );
        //  printf("  arg2: %f\n", arg2);
        while( itoken < nTokens ) {
          //    cout << "       ***  Keep on solving " << arg1 << " " << arg2 << endl;
          tokenOp2Ptr = &myExpressionList->at(iex).at(itoken++);
          arg3 = getValue( &(myExpressionList->at(iex).at(itoken++)) );
          if( tokenOp1Ptr->type == OPERATOR_PLUS_MINUS ) {
            if( tokenOp2Ptr->type == OPERATOR_MULT_DIV ) {  // Second operator is */ --> this takes precedence
              arg2 = (tokenOp2Ptr->valChar == '*') ? arg2*arg3 : arg2/arg3;     
            }
            else {  // Both operators are + -
              // Must solve arg1 & arg2 first in case next operator is */ 
              arg1 = (tokenOp1Ptr->valChar == '+') ? arg1+arg2 : arg1-arg2;
              tokenOp1Ptr = tokenOp2Ptr;
              arg2 = arg3;
            }
          }
          else {   // if( tokenOp1.type == OPERATOR_MULT_DIV ) {
            // Evaluate arg1 & arg2 --> arg1
            arg1 = (tokenOp1Ptr->valChar == '*') ? arg1*arg2 : arg1/arg2;
            tokenOp1Ptr = tokenOp2Ptr;
            arg2 = arg3;
          }
        }  // End while
        // Final evaluation of arg1 & arg2
        if( tokenOp1Ptr->type == OPERATOR_PLUS_MINUS ) {
          arg1 = (tokenOp1Ptr->valChar == '+') ? arg1+arg2 : arg1-arg2;
        }
        else {   // if( tokenOp1.type == OPERATOR_MULT_DIV ) {
          arg1 = (tokenOp1Ptr->valChar == '*') ? arg1*arg2 : arg1/arg2;
        }
      }
    }
    myValueStack[iex] = arg1;
    //        printf("  Stack result arg1: %f\n", arg1);
  }
  //  printf("  Final result: %f\n", myValueStack[size-1]);
  return myValueStack[size-1];
}

//double csEquationSolver::getValue( csToken const& token, double const* myValueStack ) {
double csEquationSolver::getValue( csToken const* token ) {
  switch( token->type ) { 
  case NUMBER:
    return token->valDouble;
  case INTERNAL_VAR:
    return myValueStack[token->valInt];
  case USER_CONSTANT:
    return myUserConstants[token->valInt];
  default:
    printf("ERROR!\n");
    return 0.0;
  }
}

//--------------------------------------------------------------------------------
// Add all checks for correctness of expression
//
void csEquationSolver::prepareExpressionList( csVector<csToken>& tokenList ) {
  csStack<csToken> tmpTokenStack;
  csVector<csToken> simpleExpression;
  csStack<csToken> functionStack;
  csToken token, tokenTmp;
  myExpressionList->clear();
  csVector<csToken> prefixExpression(3);
  tokenTmp.type = NUMBER;
  tokenTmp.valDouble = -1.0;
  prefixExpression.insertEnd(tokenTmp);
  tokenTmp.type = OPERATOR_MULT_DIV;
  tokenTmp.valChar = '*';
  prefixExpression.insertEnd(tokenTmp);
  prefixExpression.insertEnd(tokenTmp); // Dummy

  for( int i = 0; i < tokenList.size() ; i++ ) {
    token = tokenList.at(i);
    //   cout << "Token '" << token << "': ";
    if( token.type == FUNCTION ) {
      //  cout << " Push function..\n";
      functionStack.push( token );
    }
    else if( token.type == BRACKET_CLOSE || token.type == COMMA_OPERATOR ) {
      // Push enclosed simple expression into myExpressionList
      simpleExpression.clear();
      //    cout << "Push token into expression: " ;
      while( !tmpTokenStack.isEmpty() ) {
        tokenTmp = tmpTokenStack.pop();
        if( tokenTmp.type == BRACKET_OPEN ) {
          break;
        }
        else if( tokenTmp.type == COMMA_OPERATOR) {
          break;
        }
        else if( tokenTmp.type == OPERATOR_PREFIX_MIN) {
          if( simpleExpression.size() == 0 ) {
            throw( EquationException( std::string("Syntax error in equation. Unexpected prefix '-' sign.\n"), "" ) );     
          }
          if( simpleExpression.at(0).type == NUMBER ) {
            simpleExpression.at(0).valDouble *= -1.0;
          }
          else {
            //      cout << " Set prefix minus: " << simpleExpression.at(0) << endl;
            prefixExpression.at(2) = simpleExpression.at(0);
            myExpressionList->insertEnd( prefixExpression );
            tokenTmp.type = INTERNAL_VAR;
            tokenTmp.valInt = myExpressionList->size()-1;
            simpleExpression.at(0) = tokenTmp;
          }
        }
        else {
          //            cout << tokenTmp << ' ';
          simpleExpression.insert(tokenTmp,0);
        }
      }
      //          cout << endl;
      if( simpleExpression.size() > 1 ) {  // If expression has more than one item, push as separate expression into the expression 'super' list
        //      cout << "   +++++ Push the whole damn thing into expression slot #" << myExpressionList->size() << endl;
        myExpressionList->insertEnd(simpleExpression); 
        csToken t;
        t.valInt = myExpressionList->size()-1;
        t.type = INTERNAL_VAR;
        tmpTokenStack.push( t );
      }
      else if( simpleExpression.size() == 1 ) { // Expression consists of one item only --> push back into expression stack
        tmpTokenStack.push( simpleExpression.at(0) );
      }
      else {
        throw( EquationException( std::string("Empty field between tokens '") + tokenTmp.valChar + token.valChar + std::string("'.\n"), "" ) );
      }
      if( token.type == COMMA_OPERATOR ) {
        tmpTokenStack.push(token); // Keep comma for next simple expression
      }
      else if( token.type == BRACKET_CLOSE && !functionStack.isEmpty() ) {
        functionStack.pop();   // Remove open bracket
        if( functionStack.isEmpty() ) {
          throw( EquationException( std::string("Syntax error in equation.\n"), "" ) );
        }
        if( functionStack.peek().type == FUNCTION ) {
          csToken tokenFunc = functionStack.pop();
          int nArg = tokenFunc.function->nArguments;  // Assume 1 or 2 arguments in the following checks
          if( nArg <= tmpTokenStack.size() ) {
            csToken tokenArg1 = tmpTokenStack.pop(); //tmpTokenStack.at( tmpTokenStack.size()-1 );  // 
            if( !tokenArg1.isNumber() ) {
              if( tokenArg1.type == BRACKET_OPEN )
                throw( EquationException( std::string("Missing argument in math function '") + tokenFunc.function->name + std::string("'.\n"), "" ) );
              else
                throw( EquationException( std::string("Syntax error within argument of math function '") + tokenFunc.function->name + std::string("'.\n"), "" ) );
            }
            if( nArg == 2 ) {
              csToken tokenArg2 = tmpTokenStack.pop(); // tmpTokenStack.at( tmpTokenStack.size()-2 ); // tmpTokenStack.pop();
              if( !tokenArg2.isNumber() ) {
                if( tokenArg2.type == BRACKET_OPEN || tokenArg2.type == COMMA_OPERATOR )
                  throw( EquationException( std::string("Too few arguments in math function '") + tokenFunc.function->name + std::string("'.\n"), "" ) );
                else
                  throw( EquationException( std::string("Syntax error within argument of math function '") + tokenFunc.function->name + std::string("'.\n"), "" ) );
              }
              //else if( currentExpression.size() > 2 ) {
              else {
                simpleExpression.clear();
                simpleExpression.insertEnd( tokenFunc );
                simpleExpression.insertEnd( tokenArg2 );  // Actually arg2 is the first argument
                simpleExpression.insertEnd( tokenArg1 );
                myExpressionList->insertEnd( simpleExpression );
                csToken t;
                t.type = INTERNAL_VAR;
                t.valInt = myExpressionList->size()-1;
                tmpTokenStack.push( t );
              }
            }
            else {
              simpleExpression.clear();
              simpleExpression.insertEnd( tokenFunc );
              simpleExpression.insertEnd( tokenArg1 );
              myExpressionList->insertEnd( simpleExpression );
              csToken t;
              t.type = INTERNAL_VAR;
              t.valInt = myExpressionList->size()-1;
              tmpTokenStack.push( t );
            }
          }
          else {  // 
            throw( EquationException( std::string("Too few arguments in math function '") + tokenFunc.function->name + std::string("'.\n"), "" ) );
          }
        }  // END if this is a function, else nothing
      }
    }
    else {
      if( token.type == BRACKET_OPEN && !functionStack.isEmpty() ) { 
        functionStack.push(token);
      }
      tmpTokenStack.push(token);
    }
  }

  // Final bit..
  simpleExpression.clear();
  while( !tmpTokenStack.isEmpty() ) {
    tokenTmp = tmpTokenStack.pop();
    if( tokenTmp.type == BRACKET_OPEN || tokenTmp.type == COMMA_OPERATOR) {
      break;
    }
    else if( tokenTmp.type == OPERATOR_PREFIX_MIN ) {
      if( simpleExpression.size() == 0 ) {
        throw( EquationException( std::string("Syntax error in equation. Unexpected prefix '-' sign.\n"), "" ) );      
      }
      if( simpleExpression.at(0).type == NUMBER ) {
        simpleExpression.at(0).valDouble *= -1.0;
      }
      else {
        prefixExpression.at(2) = simpleExpression.at(0);
        myExpressionList->insertEnd( prefixExpression );
        tokenTmp.type = INTERNAL_VAR;
        tokenTmp.valInt = myExpressionList->size()-1;
        simpleExpression.at(0) = tokenTmp;
      }
    }
    else {
      simpleExpression.insert(tokenTmp,0);
    }
  }
  myExpressionList->insertEnd(simpleExpression);


  //------------------------------
  // Final syntax check, partially redundant. Notsure if everything is caught...
  //
  csToken* tokenPtr;
  /*
    for( int i = 0; i < myExpressionList->size(); i++ ) {
    simpleExpression = myExpressionList->at(i);
    cout << "Expression " << i << ":" << endl;
    for( int k = 0; k < simpleExpression.size(); k++ ) {
    tokenPtr = &simpleExpression.at(k);
    cout << *tokenPtr << "   ";
    }
    cout << endl;
    }
  */
  
  for( int i = 0; i < myExpressionList->size(); i++ ) {
    //    if( DEBUG_EQ_SOLVER ) cout << "Simple expression #" << i << ": ";
    simpleExpression = myExpressionList->at(i);
    if( simpleExpression.size() == 0 ) {
      throw( EquationException( std::string("Unknown error occurred in equation..\n"), "" ) );
    }
    else if( simpleExpression.at(0).type == FUNCTION ) {
      if( simpleExpression.at(0).function->nArguments != simpleExpression.size()-1 ) {
        throw( EquationException( std::string("Wrong number of arguments in math function '") + simpleExpression.at(0).function->name + std::string("'\n"), "" ) );
      }
    }
    else {
      for( int k = 0; k < simpleExpression.size(); k++ ) {
        tokenPtr = &simpleExpression.at(k);
        if( int(k/2)*2 != k ) {
          if( tokenPtr->type != OPERATOR_PLUS_MINUS && tokenPtr->type != OPERATOR_MULT_DIV ) {    
            throw( EquationException( std::string("Syntax error in equation: Missing operator +-*/ ?"), "" ) );
          }
        }
        else if( !tokenPtr->isNumber() ) {
          //    throw( ExpressionException( string("Syntax error in equation: Missing value/constant expression..?"), csToken() ) );
          throw( EquationException( std::string("Syntax error in equation: Missing value/constant expression..?"), "" ) );
        }
      }
      if( int(simpleExpression.size()/2)*2 == simpleExpression.size() ) {  // A simple expresison must have uneven number of tokens (e.q. '1+2*3')
        throw( EquationException( std::string("Syntax error in equation: Missing value/constant expression..?"), "" ) );
        //simpleExpression.at(simpleExpression.size()-1)
      }
    }
    for( int k = 0; k < simpleExpression.size(); k++ ) {
      token = simpleExpression.at(k);
      //      if( DEBUG_EQ_SOLVER ) std::cout << token << ' ';
    }
    //    if( DEBUG_EQ_SOLVER )   std::cout << endl;
  }
  if( DEBUG_EQ_SOLVER ) std::cout << "-------------------------------------\n";
  myValueStack = new double[myExpressionList->size()];
}



//-----------------------------------------------
ExpressionException::ExpressionException( std::string message ) :
  csException( message.c_str() ) {
  //  token = new csToken();
  token = NULL;
}
ExpressionException::ExpressionException( std::string message, csToken const& token_new ) :
  csException( message.c_str() ) {
  token = NULL;
  //  token = new csToken(token_new);
}
ExpressionException::ExpressionException( ExpressionException const& obj  ) : csException( obj.myMessage.c_str() ) {
}
ExpressionException::~ExpressionException() {
  if( token ) {
    //    delete token;  // Token deletion causes crash. MathFunction destructor actually causes the crash. Why?
    token = NULL;
  }
}
//-----------------------------------------------
ConstantException::ConstantException( std::string message ) :
  csException( message.c_str() ) {
}

