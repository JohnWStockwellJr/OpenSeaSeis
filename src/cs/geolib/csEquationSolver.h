/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef CS_EQUATION_SOLVER_H
#define CS_EQUATION_SOLVER_H

#include <cmath>
#include <string>
#include "csToken.h"
#include "csException.h"

namespace cseis_geolib {

// Forward declarations:
  class csMathFunction;
  class csMathConstant;
  template <typename T> class csVector;

/**
 * Parse & evaluate mathematical equations.
 *
 * This class is optimized for the case where an equation is parsed once and then solved many times
 * using different values (user constants).
 *
 * How to use:
 *  - Create new csEquationSolver object
 *  - Call prepare() method
 *  - Call prepareUserConstants() method if necessary
 *  - Call setUserConstants() method
 *  - Call solve() method
 * How it works:
 *  - prepare() parses text and prepares all internal stuff
 *  - prepareUserConstants() returns a list of all user constant names actually used in equation
 *  - setUserConstants() sets the values of the necessary user constants
 *  - solve() solves the equation and returns the result
 *
 * @author Bjorn Olofsson
 * @date 2005
 */
class csEquationSolver {
public:
  csEquationSolver();
  csEquationSolver( std::string nameSpecialFloatArray, int numValuesSpecialFloatArray );
  ~csEquationSolver();
  /**
   * Call to retrieve error message if error occurred in any of the methods (false return value)
   */
  std::string const& getErrorMessage() const;
  /**
   * Prepare equation solver
   * - Parse equation text, check user constant names etc.
   *
   * @return false if error occurred
   */
  bool prepare( std::string const& equation, std::string const* userConstantNames, int nConstants, double const* userConstants );
  bool prepare( std::string const& equation, std::string const* userConstantNames, int nConstants );
  bool prepare( std::string const& equation );
  /**
  * Return a list of all user constant names that are used in equation
  */
  void prepareUserConstants( csVector<std::string>* const constantList );

  /**
  * Set values of all user constants
  */
  void setUserConstants( double const* userConstants, int nConstants );

  /**
  * Solve equation
  * @return result
  */
  double solve();

private:
  void init();
  void prepareExpressionList( csVector<csToken>& tokenList );
  double getValue( csToken const* tokenPtr );
  /// Tokenize regular expression (lower case, no spaces)
  void tokenizeExpression( std::string expression, csVector<csToken>& tokenList );
  /// @return boolean   true if character c is valid operator
  bool isOperator(char c) const;
  /// @param expression std::string   expression containing operator to be checked for prefix notation
  /// @param operatorPos int     position in expression where operator is located
  /// @return boolean            true if operator at operatorPos is prefix operator
  bool isPrefixOperator( std::string const& expression, int operatorPos ) const;
  /// @return boolean   true if identifierName has been recognized as available constant
  bool isUserConstant( std::string const& identifierName ) const;
  int getCharType( char c ) const;
  int getUserConstantIndex( std::string const& identifierName ) const;

  csMathFunction const* getMathFunction( std::string functionName );
  csMathConstant const* getMathConstant( std::string constantName );

private:
  static bool const DEBUG_EQ_SOLVER = false;
  /// Number of user constants
  int myNumUserConstants;
  /// Values of user constants
  double* myUserConstants;    // Array of user constant values (3 rows, for index i-1, i and i+1)
  /// Identifier names for user constants
  std::string*  myUserConstantNames;
  std::string myErrorMessage;
  /// Reduced list of expressions that is solved in solve() method in order to solve the original equation.
  /// Each expression is a simple expression involving only either a math function or a succession of */+- operations
  csVector< csVector<csToken> >* myExpressionList;
  /// Stack of fixed size that is used in the solve() method.
  /// Contains solution values for each of the expressions in the above expression list
  double* myValueStack;
  int myNumValuesSpecialFloatArray;
  std::string myNameSpecialFloatArray;
  bool myIsSpecialFloatArray;
};

//------------------------------------------------------------------------------
class EquationException : public csException {
public:
  std::string equation;
  EquationException( std::string message, std::string equation_new ) :
    csException( message.c_str() ),
    equation( equation_new ) {
  }
  EquationException( EquationException const& obj  ) : csException(obj.myMessage.c_str() ) {
  }
private:
  EquationException();
};

class ExpressionException : public csException {
public:
  csToken* token;
  ExpressionException( std::string message );
  ExpressionException( std::string message, csToken const& token_new );
  ExpressionException( ExpressionException const& obj  );
  ~ExpressionException();
private:
  ExpressionException();
};
//-----------------------------------------------
class ConstantException : public csException {
public:
  ConstantException( std::string message );
};

} // END namespace

#endif

