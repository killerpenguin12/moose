//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ExpressionBuilder.h"

std::shared_ptr<ExpressionBuilder::EBTermNode> ExpressionBuilder::derivative_of;
unsigned int ExpressionBuilder::current_quad_point;

ExpressionBuilder::EBTermList
operator, (const ExpressionBuilder::EBTerm & larg, const ExpressionBuilder::EBTerm & rarg)
{
  return {larg, rarg};
}

ExpressionBuilder::EBTermList
operator, (const ExpressionBuilder::EBTerm & larg, const ExpressionBuilder::EBTermList & rargs)
{
  ExpressionBuilder::EBTermList list = {larg};
  list.insert(list.end(), rargs.begin(), rargs.end());
  return list;
}

ExpressionBuilder::EBTermList
operator, (const ExpressionBuilder::EBTermList & largs, const ExpressionBuilder::EBTerm & rarg)
{
  ExpressionBuilder::EBTermList list = largs;
  list.push_back(rarg);
  return list;
}

std::ostream &
operator<<(std::ostream & os, const ExpressionBuilder::EBTerm & term)
{
  if (term._root != NULL)
    return os << *term._root;
  else
    return os << "[NULL]";
}

std::string
ExpressionBuilder::EBSymbolNode::stringify() const
{
  return _symbol;
}

std::string
ExpressionBuilder::EBTempIDNode::stringify() const
{
  std::ostringstream s;
  s << '[' << _id << ']';
  return s.str();
}

std::string
ExpressionBuilder::EBUnaryFuncTermNode::stringify() const
{
  const char * name[] = {"sin", "cos", "tan", "abs", "log", "log2", "log10", "exp", "sinh", "cosh", "sqrt", "acos", "atan"};
  std::ostringstream s;
  s << name[_type] << '(' << *_subnode << ')';
  return s.str();
}

std::string
ExpressionBuilder::EBUnaryOpTermNode::stringify() const
{
  const char * name[] = {"-", "!"};
  std::ostringstream s;

  s << name[_type];

  if (_subnode->precedence() > precedence())
    s << '(' << *_subnode << ')';
  else
    s << *_subnode;

  return s.str();
}

std::string
ExpressionBuilder::EBBinaryFuncTermNode::stringify() const
{
  const char * name[] = {"min", "max", "atan2", "hypot", "plog"};
  std::ostringstream s;
  s << name[_type] << '(' << *_left << ',' << *_right << ')';
  return s.str();
}

std::string
ExpressionBuilder::EBBinaryOpTermNode::stringify() const
{
  const char * name[] = {"+", "-", "*", "/", "%", "^", "<", ">", "<=", ">=", "=", "!=", "&&", "||"};
  std::ostringstream s;

  if (_left->precedence() > precedence())
    s << '(' << *_left << ')';
  else
    s << *_left;

  s << name[_type];

  // these operators are left associative at equal precedence
  // (this matters for -,/,&,^ but not for + and *)
  if (_right->precedence() > precedence() ||
      (_right->precedence() == precedence() &&
       (_type == SUB || _type == DIV || _type == MOD || _type == POW)))
    s << '(' << *_right << ')';
  else
    s << *_right;

  return s.str();
}

int
ExpressionBuilder::EBBinaryOpTermNode::precedence() const
{
  switch (_type)
  {
    case ADD:
    case SUB:
      return 6;
    case MUL:
    case DIV:
    case MOD:
      return 5;
    case POW:
      return 2;
    case LESS:
    case GREATER:
    case LESSEQ:
    case GREATEREQ:
      return 8;
    case EQ:
    case NOTEQ:
      return 9;
    case OR:
    case AND:
      return 10;
  }

  mooseError("Unknown type.");
}

std::string
ExpressionBuilder::EBTernaryFuncTermNode::stringify() const
{
  const char * name[] = {"if"};
  std::ostringstream s;
  s << name[_type] << '(' << *_left << ',' << *_middle << ',' << *_right << ')';
  return s.str();
}

ExpressionBuilder::EBFunction &
ExpressionBuilder::EBFunction::operator()(const ExpressionBuilder::EBTerm & arg)
{
  this->_eval_arguments = {arg};
  return *this;
}

ExpressionBuilder::EBFunction &
ExpressionBuilder::EBFunction::operator()(const ExpressionBuilder::EBTermList & args)
{
  this->_eval_arguments = EBTermList(args);
  return *this;
}

ExpressionBuilder::EBFunction &
ExpressionBuilder::EBFunction::operator=(const ExpressionBuilder::EBTerm & term)
{
  this->_arguments = this->_eval_arguments;
  this->_term = term;
  return *this;
}

ExpressionBuilder::EBFunction &
ExpressionBuilder::EBFunction::operator=(const ExpressionBuilder::EBFunction & func)
{
  this->_arguments = this->_eval_arguments;
  this->_term = EBTerm(func);
  return *this;
}

ExpressionBuilder::EBFunction::operator ExpressionBuilder::EBTerm() const
{
  unsigned int narg = _arguments.size();
  //if (narg != _eval_arguments.size())
  //  mooseError("EBFunction is used wth a different number of arguments than it was defined with.");

  // prepare a copy of the function term to perform the substitution on
  EBTerm result(_term);

  // prepare a rule list for the substitutions
  EBSubstitutionRuleList rules;
  for (unsigned i = 0; i < narg; ++i)
    rules.push_back(new EBTermSubstitution(_arguments[i], _eval_arguments[i]));

  // perform substitution
  result.substitute(rules);

  // discard rule set
  for (unsigned i = 0; i < narg; ++i)
    delete rules[i];

  return result;
}

ExpressionBuilder::EBFunction::operator std::string() const
{
  EBTerm eval;
  eval = *this; // typecast EBFunction -> EBTerm performs a parameter substitution
  std::ostringstream s;
  s << eval;
  return s.str();
}

std::string
ExpressionBuilder::EBFunction::args()
{
  unsigned int narg = _arguments.size();
  if (narg < 1)
    return "";

  std::ostringstream s;
  s << _arguments[0];
  for (unsigned int i = 1; i < narg; ++i)
    s << ',' << _arguments[i];

  return s.str();
}

unsigned int
ExpressionBuilder::EBFunction::substitute(const EBSubstitutionRule & rule)
{
  return _term.substitute(rule);
}

unsigned int
ExpressionBuilder::EBFunction::substitute(const EBSubstitutionRuleList & rules)
{
  return _term.substitute(rules);
}

#define UNARY_FUNC_IMPLEMENT(op, OP)                                                               \
  ExpressionBuilder::EBTerm op(const ExpressionBuilder::EBTerm & term)                             \
  {                                                                                                \
    mooseAssert(term._root != NULL, "Empty term provided as argument of function " #op "()");      \
    return ExpressionBuilder::EBTerm(std::make_shared<ExpressionBuilder::EBUnaryFuncTermNode>(                   \
        term.cloneRoot(), ExpressionBuilder::EBUnaryFuncTermNode::OP));                            \
  }
UNARY_FUNC_IMPLEMENT(sin, SIN)
UNARY_FUNC_IMPLEMENT(cos, COS)
UNARY_FUNC_IMPLEMENT(tan, TAN)
UNARY_FUNC_IMPLEMENT(abs, ABS)
UNARY_FUNC_IMPLEMENT(log, LOG)
UNARY_FUNC_IMPLEMENT(log2, LOG2)
UNARY_FUNC_IMPLEMENT(log10, LOG10)
UNARY_FUNC_IMPLEMENT(exp, EXP)
UNARY_FUNC_IMPLEMENT(sinh, SINH)
UNARY_FUNC_IMPLEMENT(cosh, COSH)
UNARY_FUNC_IMPLEMENT(sqrt, SQRT)
UNARY_FUNC_IMPLEMENT(acos, ACOS)
UNARY_FUNC_IMPLEMENT(atan, ATAN)

#define BINARY_FUNC_IMPLEMENT(op, OP)                                                              \
  ExpressionBuilder::EBTerm op(const ExpressionBuilder::EBTerm & left,                             \
                               const ExpressionBuilder::EBTerm & right)                            \
  {                                                                                                \
    mooseAssert(left._root != NULL,                                                                \
                "Empty term provided as first argument of function " #op "()");                    \
    mooseAssert(right._root != NULL,                                                               \
                "Empty term provided as second argument of function " #op "()");                   \
    return ExpressionBuilder::EBTerm(std::make_shared<ExpressionBuilder::EBBinaryFuncTermNode>(                  \
        left.cloneRoot(), right.cloneRoot(), ExpressionBuilder::EBBinaryFuncTermNode::OP));        \
  }
BINARY_FUNC_IMPLEMENT(min, MIN)
BINARY_FUNC_IMPLEMENT(max, MAX)
BINARY_FUNC_IMPLEMENT(atan2, ATAN2)
BINARY_FUNC_IMPLEMENT(hypot, HYPOT)
BINARY_FUNC_IMPLEMENT(plog, PLOG)

// this is a function in ExpressionBuilder (pow) but an operator in FParser (^)
ExpressionBuilder::EBTerm
pow(const ExpressionBuilder::EBTerm & left, const ExpressionBuilder::EBTerm & right)
{
  mooseAssert(left._root != NULL, "Empty term for base of pow()");
  mooseAssert(right._root != NULL, "Empty term for exponent of pow()");
  return ExpressionBuilder::EBTerm(std::make_shared<ExpressionBuilder::EBBinaryOpTermNode>(
      left.cloneRoot(), right.cloneRoot(), ExpressionBuilder::EBBinaryOpTermNode::POW));
}

#define TERNARY_FUNC_IMPLEMENT(op, OP)                                                             \
  ExpressionBuilder::EBTerm op(const ExpressionBuilder::EBTerm & left,                             \
                               const ExpressionBuilder::EBTerm & middle,                           \
                               const ExpressionBuilder::EBTerm & right)                            \
  {                                                                                                \
    mooseAssert(left._root != NULL,                                                                \
                "Empty term provided as first argument of the ternary function " #op "()");        \
    mooseAssert(middle._root != NULL,                                                              \
                "Empty term provided as second argument of the ternary function " #op "()");       \
    mooseAssert(right._root != NULL,                                                               \
                "Empty term provided as third argument of the ternary function " #op "()");        \
    return ExpressionBuilder::EBTerm(std::make_shared<ExpressionBuilder::EBTernaryFuncTermNode>(                 \
        left.cloneRoot(),                                                                          \
        middle.cloneRoot(),                                                                        \
        right.cloneRoot(),                                                                         \
        ExpressionBuilder::EBTernaryFuncTermNode::OP));                                            \
  }

TERNARY_FUNC_IMPLEMENT(conditional, CONDITIONAL)

unsigned int
ExpressionBuilder::EBUnaryTermNode::substitute(const EBSubstitutionRuleList & rules)
{
  unsigned int nrule = rules.size();

  for (unsigned int i = 0; i < nrule; ++i)
  {
    std::shared_ptr<EBTermNode> replace = rules[i]->apply(_subnode);
    if (replace != NULL)
    {
      _subnode = replace;
      return 1;
    }
  }

  return _subnode->substitute(rules);
}

unsigned int
ExpressionBuilder::EBBinaryTermNode::substitute(const EBSubstitutionRuleList & rules)
{
  unsigned int nrule = rules.size();
  unsigned int success = 0;

  for (unsigned int i = 0; i < nrule; ++i)
  {
    std::shared_ptr<EBTermNode> replace = rules[i]->apply(_left);
    if (replace != NULL)
    {
      _left = replace;
      success = 1;
      break;
    }
  }

  if (success == 0)
    success += _left->substitute(rules);

  for (unsigned int i = 0; i < nrule; ++i)
  {
    std::shared_ptr<EBTermNode> replace = rules[i]->apply(_right);
    if (replace != NULL)
    {
      _right = replace;
      return success + 1;
    }
  }

  return success + _right->substitute(rules);
}

unsigned int
ExpressionBuilder::EBTernaryTermNode::substitute(const EBSubstitutionRuleList & rules)
{
  unsigned int nrule = rules.size();
  bool left_success = false, middle_success = false, right_success = false;
  std::shared_ptr<EBTermNode> replace;

  for (unsigned int i = 0; i < nrule; ++i)
  {
    replace = rules[i]->apply(_left);
    if (replace)
    {
      _left = replace;
      left_success = true;
      break;
    }
  }

  for (unsigned int i = 0; i < nrule; ++i)
  {
    replace = rules[i]->apply(_middle);
    if (replace)
    {
      _middle = replace;
      middle_success = true;
      break;
    }
  }

  for (unsigned int i = 0; i < nrule; ++i)
  {
    replace = rules[i]->apply(_right);
    if (replace)
    {
      _right = replace;
      right_success = true;
      break;
    }
  }

  if (!left_success)
    left_success = _left->substitute(rules);
  if (!middle_success)
    middle_success = _middle->substitute(rules);
  if (!right_success)
    right_success = _right->substitute(rules);

  return left_success + middle_success + right_success;
}

unsigned int
ExpressionBuilder::EBTerm::substitute(const EBSubstitutionRule & rule)
{
  EBSubstitutionRuleList rules(1);
  rules[0] = &rule;
  return substitute(rules);
}

unsigned int
ExpressionBuilder::EBTerm::substitute(const EBSubstitutionRuleList & rules)
{
  unsigned int nrule = rules.size();

  if (_root == NULL)
    return 0;

  for (unsigned int i = 0; i < nrule; ++i)
  {
    std::shared_ptr<EBTermNode> replace = rules[i]->apply(_root);
    if (replace != NULL)
    {
      _root = replace;
      return 1;
    }
  }

  return _root->substitute(rules);
}

ExpressionBuilder::EBTermSubstitution::EBTermSubstitution(const EBTerm & find,
                                                          const EBTerm & replace)
{
  // the expression we want to substitute (has to be a symbol node)
  const std::shared_ptr<EBSymbolNode> find_root = std::dynamic_pointer_cast<EBSymbolNode>(find.getRoot());
  if (find_root == NULL)
    mooseError("Function arguments must be pure symbols.");
  _find = find_root->stringify();

  // the term we want to substitute with
  if (replace.getRoot() != NULL)
    _replace = replace.cloneRoot();
  else
    mooseError("Trying to substitute in an empty term for ", _find);
}

std::shared_ptr<ExpressionBuilder::EBTermNode>
ExpressionBuilder::EBTermSubstitution::substitute(const EBSymbolNode & node) const
{
  if (node.stringify() == _find)
    return _replace;
  else
    return NULL;
}

std::shared_ptr<ExpressionBuilder::EBTermNode>
ExpressionBuilder::EBLogPlogSubstitution::substitute(const EBUnaryFuncTermNode & node) const
{
  if (node._type == EBUnaryFuncTermNode::LOG)
    return std::make_shared<EBBinaryFuncTermNode>(
        node.getSubnode(), _epsilon, EBBinaryFuncTermNode::PLOG);
  else
    return NULL;
}

ExpressionBuilder::EBTerm::EBTermVector
ExpressionBuilder::EBTerm::CreateEBTermVector(const std::string & var_name, unsigned int _op_num)
{
  EBTermVector vec;
  for(unsigned int i = 0; i < _op_num; ++i)
  {
    EBTerm temp((var_name + std::to_string(i)).c_str());
    vec.push_back(temp);
  }
  return vec;
}

ExpressionBuilder::EBMatrix::EBMatrix()
{
  _rowNum = 0;
  _colNum = 0;
  setSize(0, 0);
}

ExpressionBuilder::EBMatrix::EBMatrix(std::vector<std::vector <ExpressionBuilder::EBTerm> > FunctionMatrix)
{
  this->FunctionMatrix = FunctionMatrix;
  _rowNum = FunctionMatrix.size();
  _colNum = FunctionMatrix[0].size();
  checkSize();
}

ExpressionBuilder::EBMatrix::EBMatrix(unsigned int i, unsigned int j)
{
  setSize(i,j);
}

ExpressionBuilder::EBMatrix::EBMatrix(const RealTensorValue & rhs)
{
  this->setSize(3,3);
  for(unsigned int i = 0; i < 3; ++i)
    for(unsigned int j = 0; j < 3; ++j)
      (*this)[i][j] = EBTerm(rhs(i,j));
}

ExpressionBuilder::EBMatrix &
operator*(const ExpressionBuilder::EBMatrix & lhs, const ExpressionBuilder::EBMatrix & rhs)
{
  ExpressionBuilder::EBMatrix * result = new ExpressionBuilder::EBMatrix();
  result->checkMultSize(lhs,rhs);
  result->setSize(lhs.rowNum(),rhs.colNum());
  for(unsigned int i = 0; i < lhs.rowNum(); ++i)
    for(unsigned int j = 0; j < rhs.rowNum(); ++j)
      for(unsigned int k = 0; k < rhs.colNum(); ++ k)
        (*result)[i][k] = (*result)[i][k] + lhs[i][j] * rhs[j][k];
  return *result;
}

ExpressionBuilder::EBMatrix &
operator*(const Real & lhs, const ExpressionBuilder::EBMatrix & rhs)
{
  ExpressionBuilder::EBMatrix * result = new ExpressionBuilder::EBMatrix();
  result->setSize(rhs.rowNum(),rhs.colNum());
  for(unsigned int i = 0; i < rhs.rowNum(); ++i)
    for(unsigned int j = 0; j < rhs.colNum(); ++j)
      (*result)[i][j] = lhs * rhs[i][j];
  return *result;
}

ExpressionBuilder::EBMatrix &
operator*(const ExpressionBuilder::EBMatrix & lhs, const Real & rhs)
{
  return rhs * lhs;
}

ExpressionBuilder::EBMatrix &
operator*(const ExpressionBuilder::EBTerm & lhs, const ExpressionBuilder::EBMatrix & rhs)
{
  ExpressionBuilder::EBMatrix * result = new ExpressionBuilder::EBMatrix();
  result->setSize(rhs.rowNum(),rhs.colNum());
  for(unsigned int i = 0; i < rhs.rowNum(); ++i)
    for(unsigned int j = 0; j < rhs.colNum(); ++j)
      (*result)[i][j] = lhs * rhs[i][j];
  return *result;
}

ExpressionBuilder::EBMatrix &
operator*(const ExpressionBuilder::EBMatrix & lhs, const ExpressionBuilder::EBTerm & rhs)
{
  return rhs * lhs;
}

ExpressionBuilder::EBMatrix &
operator/(const ExpressionBuilder::EBMatrix & lhs, const ExpressionBuilder::EBTerm & rhs)
{
  ExpressionBuilder::EBMatrix * result = new ExpressionBuilder::EBMatrix();
  result->setSize(lhs.rowNum(),lhs.colNum());
  for(unsigned int i = 0; i < lhs.rowNum(); ++i)
    for(unsigned int j = 0; j < lhs.colNum(); ++j)
      (*result)[i][j] = lhs[i][j] / rhs;
  return *result;
}

ExpressionBuilder::EBMatrix &
operator/(const ExpressionBuilder::EBMatrix & lhs, const Real & rhs)
{
  ExpressionBuilder::EBMatrix * result = new ExpressionBuilder::EBMatrix();
  result->setSize(lhs.rowNum(),lhs.colNum());
  for(unsigned int i = 0; i < lhs.rowNum(); ++i)
    for(unsigned int j = 0; j < lhs.colNum(); ++j)
      (*result)[i][j] = (1/rhs) * lhs[i][j];
  return *result;
}

ExpressionBuilder::EBMatrix &
operator+(const ExpressionBuilder::EBMatrix & lhs, const ExpressionBuilder::EBMatrix & rhs)
{
  ExpressionBuilder::EBMatrix * result = new ExpressionBuilder::EBMatrix();
  result->checkAddSize(rhs,lhs);
  result->setSize(rhs.rowNum(),rhs.colNum());
  for(unsigned int i = 0; i < rhs.rowNum(); ++i)
    for(unsigned int j = 0; j < rhs.colNum(); ++j)
      (*result)[i][j] = rhs[i][j] + lhs[i][j];
  return *result;
}

ExpressionBuilder::EBMatrix &
operator-(const ExpressionBuilder::EBMatrix & lhs, const ExpressionBuilder::EBMatrix & rhs)
{
  return lhs + (-1 * rhs);
}

ExpressionBuilder::EBMatrix::operator ExpressionBuilder::EBQuaternion() const
{
  Real _epsilon = 0.0001;
  EBQuaternion newQuat;
  ExpressionBuilder::EBTerm temp = (*this)[0][0] + (*this)[1][1] + (*this)[2][2];

  newQuat[0] = conditional(temp > -(1.0 - _epsilon), sqrt(1.0 + temp) / 2.0, EBTerm(0.0));
  newQuat[1] = conditional(temp > -(1.0 - _epsilon), ((*this)[1][2] - (*this)[2][1]) / (4.0 * newQuat[0]), ExpressionBuilder::EBTerm(0.0));
  newQuat[2] = conditional(temp > -(1.0 - _epsilon), ((*this)[2][0] - (*this)[0][2]) / (4.0 * newQuat[0]), ExpressionBuilder::EBTerm(0.0));
  newQuat[3] = conditional(temp > -(1.0 - _epsilon), ((*this)[0][1] - (*this)[1][0]) / (4.0 * newQuat[0]), sqrt(-((*this)[0][0] + (*this)[1][1]) / 2.0));

  newQuat[1] = conditional(abs(newQuat[3] > _epsilon),(*this)[0][2] / (2.0 * newQuat[3]),sqrt(((*this)[0][0] + 1.0) / 2.0));
  newQuat[2] = conditional(abs(newQuat[3] > _epsilon),(*this)[1][2] / (2.0 * newQuat[3]),conditional(newQuat[1] == 0, EBTerm(1.0), (*this)[1][0] / (2.0 * newQuat[1])));
  newQuat[3] = conditional(abs(newQuat[3] > _epsilon),newQuat[3],EBTerm(0));

  newQuat[1] = -newQuat[1];
  newQuat[2] = -newQuat[2];
  newQuat[3] = -newQuat[3];

  return newQuat;
}

ExpressionBuilder::EBMatrix::operator std::vector<std::vector<std::string> >() const
{
  std::vector<std::vector<std::string> > sVec;
  for(unsigned int i = 0; i < _rowNum; ++i)
  {
    std::vector<std::string> tempVec;
    for(unsigned int j = 0; j < _colNum; ++j)
      tempVec.push_back((*this)[i][j]);
    sVec.push_back(tempVec);
  }
  return sVec;
}

ExpressionBuilder::EBMatrix &
ExpressionBuilder::EBMatrix::operator+=(const ExpressionBuilder::EBMatrix & rhs)
{
  return *this + rhs;
}

ExpressionBuilder::EBMatrix &
ExpressionBuilder::EBMatrix::operator-()
{
  return (-1) * (*this);
}

std::vector <ExpressionBuilder::EBTerm> &
ExpressionBuilder::EBMatrix::operator[](unsigned int i)
{
  if(i > _rowNum)
  {} //MooseError

  return FunctionMatrix[i];
}

const std::vector <ExpressionBuilder::EBTerm> &
ExpressionBuilder::EBMatrix::operator[](unsigned int i) const
{
  if(i > _rowNum)
  {} //MooseError

  return FunctionMatrix[i];
}

ExpressionBuilder::EBMatrix
ExpressionBuilder::EBMatrix::transpose()
{
  for(unsigned int i = 0; i < _rowNum; ++i)
    for(unsigned int j = i + 1; j < _colNum; ++j)
    {
      ExpressionBuilder::EBTerm temp = (*this)[i][j];
      (*this)[i][j] = (*this)[j][i];
      (*this)[j][i] = temp;
    }
  return *this;
}

void
ExpressionBuilder::EBMatrix::checkSize()
{
  for(unsigned int i = 0; i < _rowNum; ++i)
    if(FunctionMatrix[i].size() != _colNum)
    {} //MooseError
}

unsigned int
ExpressionBuilder::EBMatrix::rowNum() const
{
  return _rowNum;
}

unsigned int
ExpressionBuilder::EBMatrix::colNum() const
{
  return _colNum;
}

void
ExpressionBuilder::EBMatrix::setSize(const unsigned int i, const unsigned int j)
{
  _rowNum = i;
  _colNum = j;
  FunctionMatrix.resize(i);
  for(unsigned int k = 0; k < i; k++)
    for(unsigned int l = 0; l < j; ++l)
      FunctionMatrix[k].push_back(EBTerm(0));
}

ExpressionBuilder::EBMatrix
ExpressionBuilder::EBMatrix::rotVec1ToVec2(ExpressionBuilder::EBVector & Vec1, ExpressionBuilder::EBVector & Vec2)
{
  EBMatrix rot1_to_z = rotVecToZ(Vec1);
  EBMatrix rot2_to_z = rotVecToZ(Vec2);
  return (rot2_to_z.transpose()) * rot1_to_z;
}

ExpressionBuilder::EBMatrix
ExpressionBuilder::EBMatrix::rotVecToZ(ExpressionBuilder::EBVector & Vec)
{
  Vec = Vec / Vec.ExpressionBuilder::EBVector::norm();

  ExpressionBuilder::EBVector v1;
  ExpressionBuilder::EBVector w;

  w[0] = abs(Vec[0]);
  w[1] = abs(Vec[1]);
  w[2] = abs(Vec[2]);

  v1[0] = conditional((w[2] >= w[1] && w[1] >= w[0]) || (w[1] >= w[2] && w[2] >= w[0]), EBTerm(1.0),EBTerm(0.0));
  v1[1] = conditional((w[2] >= w[0] && w[0] >= w[1]) || (w[0] >= w[2] && w[2] >= w[1]), EBTerm(1.0),EBTerm(0.0));
  v1[2] = conditional((v1[1] == 0.0 && v1[1] == 0.0),1.0,0.0);

  v1 = v1 - (v1 * Vec) * Vec;
  v1 = v1 / v1.ExpressionBuilder::EBVector::norm();

  ExpressionBuilder::EBVector v0;
  v0[0] = v1[1] * Vec[2] - v1[2] * Vec[1];
  v0[1] = v1[2] * Vec[0] - v1[0] * Vec[2];
  v0[2] = v1[0] * Vec[1] - v1[1] * Vec[0];

  ExpressionBuilder::EBMatrix newMat(3,3);

  newMat[0][0] = v0[0];
  newMat[0][1] = v0[1];
  newMat[0][2] = v0[2];
  newMat[1][0] = v1[0];
  newMat[1][1] = v1[1];
  newMat[1][2] = v1[2];
  newMat[2][0] = Vec[0];
  newMat[2][1] = Vec[1];
  newMat[2][2] = Vec[2];

  return newMat;
}

void
ExpressionBuilder::EBMatrix::checkMultSize(const ExpressionBuilder::EBMatrix & lhs, const ExpressionBuilder::EBMatrix & rhs)
{
  if(lhs.colNum() != rhs.rowNum())
  {} //MooseError
}

void
ExpressionBuilder::EBMatrix::checkAddSize(const ExpressionBuilder::EBMatrix & lhs, const ExpressionBuilder::EBMatrix & rhs)
{
  if(lhs.colNum() != rhs.colNum())
  {} //MooseError
  if(lhs.rowNum() != rhs.rowNum())
  {} //MooseError
}

ExpressionBuilder::EBVector::EBVector()
{
  for(int i = 0; i < 3; ++i)
  {
    FunctionVector.push_back(EBTerm(0));
  }
}

ExpressionBuilder::EBVector::EBVector(std::vector<EBTerm> FunctionVector)
{
  checkSize(FunctionVector);
  this->FunctionVector = FunctionVector;
}

ExpressionBuilder::EBVector::EBVector(std::vector<std::string> FunctionNameVector)
{
  checkSize(FunctionNameVector);
  for(int i = 0; i < 3; ++i)
    this->FunctionVector.push_back(EBTerm(FunctionNameVector[i].c_str()));
}

ExpressionBuilder::EBVector &
ExpressionBuilder::EBVector::operator=(const RealVectorValue & rhs)
{
  EBVector * result = new EBVector;
  for(unsigned int i = 0; i < 3; ++i)
  {
    (*result)[i] = EBTerm(rhs(i));
  }
  return (*result);
}

ExpressionBuilder::EBVector::EBVectorVector
ExpressionBuilder::EBVector::CreateEBVectorVector(const std::string & var_name, unsigned int _op_num)
{
  EBVectorVector vec;
  for(unsigned int i = 0; i < _op_num; ++i)
  {
    EBVector temp;
    EBTerm tempx((var_name + std::to_string(i) + "x").c_str());
    EBTerm tempy((var_name + std::to_string(i) + "y").c_str());
    EBTerm tempz((var_name + std::to_string(i) = "z").c_str());
    temp.push_back(tempx);
    temp.push_back(tempy);
    temp.push_back(tempz);
    vec.push_back(temp);
  }
  return vec;
}

ExpressionBuilder::EBVector::operator std::vector<std::string>() const
{
  std::vector<std::string> sVec;
  for(unsigned int i = 0; i < 3; ++i)
    sVec.push_back((*this)[i]);
  return sVec;
}

ExpressionBuilder::EBTerm &
operator*(const ExpressionBuilder::EBVector & lhs, const ExpressionBuilder::EBVector & rhs)
{
  ExpressionBuilder::EBTerm * result = new ExpressionBuilder::EBTerm(0);
  for(unsigned int i = 0; i < 3; ++i)
    *result = *result + lhs[i] * rhs[i];
  return *result;
}

ExpressionBuilder::EBVector &
operator*(const ExpressionBuilder::EBVector & lhs, const Real & rhs)
{
  ExpressionBuilder::EBVector * result = new ExpressionBuilder::EBVector;
  for(unsigned int i = 0; i < 3; ++i)
    (*result)[i] = rhs * lhs[i];
  return *result;
}

ExpressionBuilder::EBVector &
operator*(const Real & lhs, const ExpressionBuilder::EBVector & rhs)
{
  return rhs * lhs;
}

ExpressionBuilder::EBVector &
operator*(const ExpressionBuilder::EBTerm & lhs, const ExpressionBuilder::EBVector & rhs)
{
  ExpressionBuilder::EBVector * result = new ExpressionBuilder::EBVector;
  for(unsigned int i = 0; i < 3; ++i)
    (*result)[i] = rhs[i] * lhs;
  return *result;
}

ExpressionBuilder::EBVector &
operator*(const ExpressionBuilder::EBVector & lhs, const ExpressionBuilder::EBTerm & rhs)
{
  return rhs * lhs;
}

ExpressionBuilder::EBVector &
operator*(const ExpressionBuilder::EBVector & lhs, const ExpressionBuilder::EBMatrix & rhs) // We assume the vector is 1 x 3 here
{
  if(rhs.rowNum() != 3)
  {} // MooseError
  ExpressionBuilder::EBVector * result = new ExpressionBuilder::EBVector;
  for(unsigned int i = 0; i < 3; ++i)
    for(unsigned int j =0; j < 3; ++j)
      (*result)[i] = (*result)[i] + lhs[j] * rhs[j][i];
  return *result;
}

ExpressionBuilder::EBVector &
operator*(const ExpressionBuilder::EBMatrix & lhs, const ExpressionBuilder::EBVector & rhs) // We assume the vector is 3 x 1 here
{
  if(lhs.colNum() != 3)
  {} // MooseError
  ExpressionBuilder::EBVector * result = new ExpressionBuilder::EBVector;
  for(unsigned int i = 0; i < 3; ++i)
    for(unsigned int j =0; j < 3; ++j)
      (*result)[i] = (*result)[i] + lhs[i][j] * rhs[j];
  return *result;
}

ExpressionBuilder::EBVector &
operator/(const ExpressionBuilder::EBVector & lhs, const Real & rhs)
{
  ExpressionBuilder::EBVector * result = new ExpressionBuilder::EBVector;
  for(unsigned int i = 0; i < 3; ++i)
  {
    (*result)[i] = (1/rhs) * lhs[i];
  }
  return *result;
}

ExpressionBuilder::EBVector &
operator/(const ExpressionBuilder::EBVector & lhs, const ExpressionBuilder::EBTerm & rhs)
{
  ExpressionBuilder::EBVector * result = new ExpressionBuilder::EBVector;
  for(unsigned int i = 0; i < 3; ++i)
    (*result)[i] = lhs[i] / rhs;
  return (*result);
}

ExpressionBuilder::EBVector &
operator+(const ExpressionBuilder::EBVector & lhs, const ExpressionBuilder::EBVector & rhs)
{
  ExpressionBuilder::EBVector * result = new ExpressionBuilder::EBVector;
  for(unsigned int i = 0; i < 3; ++i)
    (*result)[i] = lhs[i] + rhs[i];
  return *result;
}

ExpressionBuilder::EBVector &
operator-(const ExpressionBuilder::EBVector & lhs, const ExpressionBuilder::EBVector & rhs)
{
  return lhs + ((-1) * rhs);
}

ExpressionBuilder::EBVector &
ExpressionBuilder::EBVector::operator+=(const ExpressionBuilder::EBVector & rhs)
{
  return rhs + *this;
}

ExpressionBuilder::EBVector &
ExpressionBuilder::EBVector::operator-()
{
  return (-1) * (*this);
}

ExpressionBuilder::EBTerm &
ExpressionBuilder::EBVector::operator[](unsigned int i)
{
  return FunctionVector[i];
}

const ExpressionBuilder::EBTerm &
ExpressionBuilder::EBVector::operator[](unsigned int i) const
{
  return FunctionVector[i];
}

ExpressionBuilder::EBVector
ExpressionBuilder::EBVector::cross(const ExpressionBuilder::EBVector & lhs, const ExpressionBuilder::EBVector & rhs)
{
  ExpressionBuilder::EBVector * result = new ExpressionBuilder::EBVector;
  (*result)[0] = lhs[1] * rhs[2] - lhs[2] * rhs[1];
  (*result)[1] = lhs[0] * rhs[2] - lhs[2] * rhs[0];
  (*result)[2] = lhs[0] * rhs[1] - lhs[1] * rhs[0];
  return *result;
}

ExpressionBuilder::EBTerm
ExpressionBuilder::EBVector::norm()
{
  return sqrt(*this * *this);
}

void
ExpressionBuilder::EBVector::push_back(EBTerm term)
{
  FunctionVector.push_back(term);
  if(FunctionVector.size() > 3)
  {} //MooseError
}

void
ExpressionBuilder::EBVector::checkSize(std::vector<EBTerm> FunctionVector)
{
  if(FunctionVector.size() != 3)
  {} // MooseError
}

void
ExpressionBuilder::EBVector::checkSize(std::vector<std::string> FunctionVector)
{
  if(FunctionVector.size() != 3)
  {} // MooseError
}

ExpressionBuilder::EBQuaternion::EBQuaternion(std::vector<ExpressionBuilder::EBTerm> FunctionQuat)
{
  this->FunctionQuat = FunctionQuat;
}

ExpressionBuilder::EBQuaternion &
operator*(const ExpressionBuilder::EBQuaternion & lhs, const ExpressionBuilder::EBQuaternion & rhs)
{
  ExpressionBuilder::EBQuaternion * result = new ExpressionBuilder::EBQuaternion;
  (*result)[0] = lhs[0] * rhs[0] - lhs[1] * rhs[1] - lhs[2] * rhs[2] - lhs[3] * rhs[3];
  (*result)[1] = lhs[0] * rhs[1] + lhs[1] * rhs[0] + lhs[2] * rhs[3] - lhs[3] * rhs[2];
  (*result)[2] = lhs[0] * rhs[2] - lhs[1] * rhs[3] + lhs[2] * rhs[0] + lhs[3] * rhs[1];
  (*result)[3] = lhs[0] * rhs[3] + lhs[1] * rhs[2] - lhs[2] * rhs[1] + lhs[3] * rhs[0];
  return *result;
}

ExpressionBuilder::EBQuaternion &
operator*(const ExpressionBuilder::EBQuaternion & lhs, const Real & rhs)
{
  ExpressionBuilder::EBQuaternion * result = new ExpressionBuilder::EBQuaternion;
  for(unsigned int i = 0; i < 4; ++i)
  {
    (*result)[i] = lhs[i] * rhs;
  }
  return (*result);
}

ExpressionBuilder::EBQuaternion &
operator*(const Real & lhs, const ExpressionBuilder::EBQuaternion & rhs)
{
  return rhs * lhs;
}

ExpressionBuilder::EBQuaternion &
operator*(const ExpressionBuilder::EBQuaternion & lhs, const ExpressionBuilder::EBTerm & rhs)
{
  ExpressionBuilder::EBQuaternion * result = new ExpressionBuilder::EBQuaternion;
  for(unsigned int i = 0; i < 4; ++i)
  {
    (*result)[i] = lhs[i] * rhs;
  }
  return (*result);
}

ExpressionBuilder::EBQuaternion &
operator*(const ExpressionBuilder::EBTerm & lhs, const ExpressionBuilder::EBQuaternion & rhs)
{
  return rhs * lhs;
}

ExpressionBuilder::EBQuaternion &
operator/(const ExpressionBuilder::EBQuaternion & lhs, const Real & rhs)
{
  return (1/rhs) * lhs;
}

ExpressionBuilder::EBQuaternion &
operator/(const ExpressionBuilder::EBQuaternion & lhs, const ExpressionBuilder::EBTerm & rhs)
{
  return (1/rhs) * lhs;
}

ExpressionBuilder::EBQuaternion &
operator+(const ExpressionBuilder::EBQuaternion & lhs, const ExpressionBuilder::EBQuaternion & rhs)
{
  ExpressionBuilder::EBQuaternion * result = new ExpressionBuilder::EBQuaternion;
  for(unsigned int i = 0; i < 4; ++i)
  {
    (*result)[i] = lhs[i] + rhs[i];
  }
  return (*result);
}

ExpressionBuilder::EBQuaternion &
operator-(const ExpressionBuilder::EBQuaternion & lhs, const ExpressionBuilder::EBQuaternion & rhs)
{
  return lhs + ((-1) * rhs);
}

ExpressionBuilder::EBQuaternion &
ExpressionBuilder::EBQuaternion::operator-(const ExpressionBuilder::EBQuaternion & rhs)
{
  return (-1) * rhs;
}

ExpressionBuilder::EBTerm &
ExpressionBuilder::EBQuaternion::operator[](unsigned int i)
{
  return FunctionQuat[i];
}

const ExpressionBuilder::EBTerm &
ExpressionBuilder::EBQuaternion::operator[](unsigned int i) const
{
  return FunctionQuat[i];
}

ExpressionBuilder::EBQuaternion::operator ExpressionBuilder::EBMatrix() const
{
  EBMatrix newMat(3,3);
  newMat[0][0] = (*this)[0] * (*this)[0] + (*this)[1] * (*this)[1] - (*this)[2] * (*this)[2] - (*this)[3] * (*this)[3];
  newMat[0][1] = 2 * (*this)[1] * (*this)[2] - 2 * (*this)[0] * (*this)[3];
  newMat[0][2] = 2 * (*this)[1] * (*this)[3] + 2 * (*this)[0] * (*this)[2];
  newMat[1][0] = 2 * (*this)[1] * (*this)[2] + 2 * (*this)[0] * (*this)[3];
  newMat[1][1] = (*this)[0] * (*this)[0] - (*this)[1] * (*this)[1] + (*this)[2] * (*this)[2] - (*this)[3] * (*this)[3];
  newMat[1][2] = 2 * (*this)[2] * (*this)[3] - 2 * (*this)[0] * (*this)[1];
  newMat[2][0] = 2 * (*this)[1] * (*this)[3] - 2 * (*this)[0] * (*this)[2];
  newMat[2][1] = 2 * (*this)[2] * (*this)[3] + 2 * (*this)[0] * (*this)[1];
  newMat[2][2] = (*this)[0] * (*this)[0] - (*this)[1] * (*this)[1] - (*this)[2] * (*this)[2] + (*this)[3] * (*this)[3];

  const EBTerm d = (*this)[0] * (*this)[0] + (*this)[1] * (*this)[1] + (*this)[2] * (*this)[2] + (*this)[3] * (*this)[3];

  return newMat / d;
}

ExpressionBuilder::EBTerm
ExpressionBuilder::EBQuaternion::norm()
{
  ExpressionBuilder::EBTerm temp;
  temp = (*this)[0] * (*this)[0];
  for(unsigned int i = 1; i < 4; ++i)
  {
    temp += (*this)[i] * (*this)[i];
  }
  return sqrt(temp);
}

void
ExpressionBuilder::EBQuaternion::checkSize(std::vector<EBTerm> FunctionQuat)
{
  if(FunctionQuat.size() != 4)
  {} //MooseError
}

/******************
 *
 * Simplification and Derivative Rules
 *
 ******************/

std::shared_ptr<ExpressionBuilder::EBTermNode>
ExpressionBuilder::EBUnaryFuncTermNode::simplify()
{
  if(!isSimplified)
  {
    std::shared_ptr<EBTermNode> new_sub = _subnode->simplify();
    if(new_sub == NULL)
    {
      _subnode = std::make_shared<EBNumberNode<int> >(0);
      if(_type == SIN || _type == TAN || _type == ABS || _type == SINH || _type == SQRT || _type == ATAN)
      {
        simplified = NULL;
        isSimplified = true;
        return simplified;
      }
    }
    else if(_subnode != new_sub)
      _subnode = new_sub;
    simplified = this->shared_from_this();
    isSimplified = true;
  }
  return simplified;
}

std::shared_ptr<ExpressionBuilder::EBTermNode>
ExpressionBuilder::EBBinaryOpTermNode::simplify()
{
  if(!isSimplified)
  {
    std::shared_ptr<EBTermNode> new_left = _left->simplify();
    std::shared_ptr<EBTermNode> new_right = _right->simplify();
    if(new_left == NULL)
    {
      _left = std::make_shared<EBNumberNode<int> >(0);
      if(_type == MUL || _type == DIV)
      {
        simplified = NULL;
        isSimplified = true;
        return simplified;
      }
      if(new_right == NULL)
      {
        _right = std::make_shared<EBNumberNode<int> >(0);
        if(_type == ADD || _type == SUB)
        {
          simplified = NULL;
          isSimplified = true;
          return simplified;
        }
      }
      if(_type == ADD)
      {
        simplified = _right;
        isSimplified = true;
        return simplified;
      }
    }
    else if(_left != new_left)
      _left = new_left;
    if(new_right == NULL)
    {
      _right = std::make_shared<EBNumberNode<int> >(0);
      if(_type == MUL)
      {
        simplified = NULL;
        isSimplified = true;
        return simplified;
      }
      if(_type == ADD || _type == SUB)
      {
        simplified = _left;
        isSimplified = true;
        return simplified;
      }
    }
    else if(_right != new_right)
      _right = new_right;
    simplified = this->shared_from_this();
    isSimplified = true;
  }
  return simplified;
}

std::shared_ptr<ExpressionBuilder::EBTermNode>
ExpressionBuilder::EBUnaryFuncTermNode::takeDerivative()
{
  if(derivative == NULL)
  {
    std::shared_ptr<EBTermNode> temp;
    switch (_type)
    {
      case SIN:
      {
        temp = std::make_shared<EBUnaryFuncTermNode>(_subnode, NodeType::COS);
        break;
      }
      case COS:
      {
        std::shared_ptr<EBUnaryFuncTermNode> temp2 = std::make_shared<EBUnaryFuncTermNode>(_subnode, NodeType::SIN);
        temp = std::make_shared<EBUnaryOpTermNode>(temp2, EBUnaryOpTermNode::NodeType::NEG);
        break;
      }
      case TAN:
      {
        std::shared_ptr<EBUnaryFuncTermNode> temp3 = std::make_shared<EBUnaryFuncTermNode>(_subnode, NodeType::COS);
        std::shared_ptr<EBBinaryOpTermNode> temp2 = std::make_shared<EBBinaryOpTermNode>(temp3, temp3, EBBinaryOpTermNode::NodeType::MUL);
        temp = std::make_shared<EBBinaryOpTermNode>(std::make_shared<EBNumberNode<Real> >(1.0), temp2, EBBinaryOpTermNode::NodeType::DIV);
        break;
      }
      case ABS:
      {
        temp = std::make_shared<EBBinaryOpTermNode>(_subnode, this->shared_from_this(), EBBinaryOpTermNode::NodeType::DIV);
        break;
      }
      case LOG:
      {
        temp = std::make_shared<EBBinaryOpTermNode>(std::make_shared<EBNumberNode<Real> >(1.0), _subnode, EBBinaryOpTermNode::NodeType::DIV);
        break;
      }
      case LOG2:
      {
        std::shared_ptr<EBBinaryOpTermNode> temp2 = std::make_shared<EBBinaryOpTermNode>(std::make_shared<EBNumberNode<Real> >(1.0), _subnode, EBBinaryOpTermNode::NodeType::DIV);
        temp = std::make_shared<EBBinaryOpTermNode>(std::make_shared<EBNumberNode<Real> >(1.0 / log(2.0)), temp2, EBBinaryOpTermNode::NodeType::MUL);
        break;
      }
      case LOG10:
      {
        std::shared_ptr<EBBinaryOpTermNode> temp2 = std::make_shared<EBBinaryOpTermNode>(std::make_shared<EBNumberNode<Real> >(1.0), _subnode, EBBinaryOpTermNode::NodeType::DIV);
        temp = std::make_shared<EBBinaryOpTermNode>(std::make_shared<EBNumberNode<Real> >(1.0 / log(10.0)), temp2, EBBinaryOpTermNode::NodeType::MUL);
        break;
      }
      case EXP:
      {
        temp = this->shared_from_this();
        break;
      }
      case SINH:
      {
        temp = std::make_shared<EBUnaryFuncTermNode>(_subnode, NodeType::COSH);
        break;
      }
      case COSH:
      {
        temp = std::make_shared<EBUnaryFuncTermNode>(_subnode, NodeType::SINH);
        break;
      }
      case SQRT:
      {
        std::shared_ptr<EBBinaryOpTermNode> temp2 = std::make_shared<EBBinaryOpTermNode>(this->shared_from_this(), std::make_shared<EBNumberNode<Real> >(2.0), EBBinaryOpTermNode::NodeType::MUL);
        temp = std::make_shared<EBBinaryOpTermNode>(std::make_shared<EBNumberNode<Real> >(1.0), temp2, EBBinaryOpTermNode::NodeType::DIV);
        break;
      }
      case ACOS:
      {
        std::shared_ptr<EBBinaryOpTermNode> temp4  = std::make_shared<EBBinaryOpTermNode>(_subnode, _subnode, EBBinaryOpTermNode::NodeType::MUL);
        std::shared_ptr<EBBinaryOpTermNode> temp3 = std::make_shared<EBBinaryOpTermNode>(std::make_shared<EBNumberNode<Real> >(1.0), temp4, EBBinaryOpTermNode::NodeType::SUB);
        std::shared_ptr<EBBinaryOpTermNode> temp2 = std::make_shared<EBBinaryOpTermNode>(std::make_shared<EBNumberNode<Real> >(1.0), temp3, EBBinaryOpTermNode::NodeType::DIV);
        temp = std::make_shared<EBUnaryOpTermNode>(temp2, EBUnaryOpTermNode::NodeType::NEG);
        break;
      }
      case ATAN:
      {
        std::shared_ptr<EBBinaryOpTermNode> temp3  = std::make_shared<EBBinaryOpTermNode>(_subnode, _subnode, EBBinaryOpTermNode::NodeType::MUL);
        std::shared_ptr<EBBinaryOpTermNode> temp2 = std::make_shared<EBBinaryOpTermNode>(std::make_shared<EBNumberNode<Real> >(1.0), temp3, EBBinaryOpTermNode::NodeType::ADD);
        temp = std::make_shared<EBBinaryOpTermNode>(std::make_shared<EBNumberNode<Real> >(1.0), temp2, EBBinaryOpTermNode::NodeType::DIV);
        break;
      }
    }
    derivative = std::make_shared<EBBinaryOpTermNode>(temp, _subnode->takeDerivative(), EBBinaryOpTermNode::NodeType::MUL);
  }
  return derivative;
}

std::shared_ptr<ExpressionBuilder::EBTermNode>
ExpressionBuilder::EBBinaryOpTermNode::takeDerivative()
{
  if(derivative == NULL)
  {
    switch (_type)
    {
      case ADD:
      {
        derivative = std::make_shared<EBBinaryOpTermNode>(_left->takeDerivative(), _right->takeDerivative(), NodeType::ADD);
        break;
      }
      case SUB:
      {
        derivative = std::make_shared<EBBinaryOpTermNode>(_left->takeDerivative(), _right->takeDerivative(), NodeType::SUB);
        break;
      }
      case MUL:
      {
        std::shared_ptr<EBBinaryOpTermNode>temp = std::make_shared<EBBinaryOpTermNode>(_left->takeDerivative(), _right, NodeType::MUL);
        std::shared_ptr<EBBinaryOpTermNode>temp2 = std::make_shared<EBBinaryOpTermNode>(_left, _right->takeDerivative(), NodeType::MUL);
        derivative = std::make_shared<EBBinaryOpTermNode>(temp, temp2, NodeType::ADD);
        break;
      }
      case DIV:
      {
        std::shared_ptr<EBBinaryOpTermNode>temp = std::make_shared<EBBinaryOpTermNode>(_right, _left->takeDerivative(), NodeType::MUL);
        std::shared_ptr<EBBinaryOpTermNode>temp2 = std::make_shared<EBBinaryOpTermNode>(_left, _right->takeDerivative(), NodeType::MUL);
        std::shared_ptr<EBBinaryOpTermNode>temp3 = std::make_shared<EBBinaryOpTermNode>(temp, temp2, NodeType::SUB);
        std::shared_ptr<EBBinaryOpTermNode>temp4 = std::make_shared<EBBinaryOpTermNode>(_right, _right, NodeType::MUL);
        derivative = std::make_shared<EBBinaryOpTermNode>(temp3, temp4, DIV);
        break;
      }
      case MOD:
      {
        derivative = _left->takeDerivative();
        break;
      }
      case POW:
      {
        std::shared_ptr<EBBinaryOpTermNode> temp = std::make_shared<EBBinaryOpTermNode>(_right, std::make_shared<EBNumberNode<Real> >(1.0), NodeType::SUB);
        std::shared_ptr<EBBinaryOpTermNode>temp2 = std::make_shared<EBBinaryOpTermNode>(_left, temp, NodeType::POW);
        derivative = std::make_shared<EBBinaryOpTermNode>(temp2, _right, MUL);
        break;
      }
      case LESS:
      case GREATER:
      case LESSEQ:
      case GREATEREQ:
      case EQ:
      case NOTEQ:
      case AND:
      case OR:
      {
        mooseError("Logical Operator not in logical statement");
        break;
      }
    }
  }
  return derivative;
}

std::shared_ptr<ExpressionBuilder::EBTermNode>
ExpressionBuilder::EBBinaryFuncTermNode::takeDerivative()
{
  if(derivative == NULL)
  {
    std::shared_ptr<EBTermNode> temp;
    std::shared_ptr<EBTermNode> temp2;
    std::shared_ptr<EBTermNode> temp3;
    std::shared_ptr<EBTermNode> temp4;
    std::shared_ptr<EBTermNode> temp5;
    std::shared_ptr<EBTermNode> temp6;
    std::shared_ptr<EBTermNode> temp7;
    std::shared_ptr<EBTermNode> temp8;
    switch (_type)
    {
      case MIN:
      {
        temp = std::make_shared<EBBinaryOpTermNode>(_left, _right, EBBinaryOpTermNode::NodeType::LESS);
        temp2 = std::make_shared<EBTernaryFuncTermNode>(temp, _left, _right, EBTernaryFuncTermNode::NodeType::CONDITIONAL);
        derivative = temp2->takeDerivative();
        break;
      }
      case MAX:
      {
        temp = std::make_shared<EBBinaryOpTermNode>(_left, _right, EBBinaryOpTermNode::NodeType::GREATER);
        temp2 = std::make_shared<EBTernaryFuncTermNode>(temp, _left, _right, EBTernaryFuncTermNode::NodeType::CONDITIONAL);
        derivative = temp2->takeDerivative();
        break;
      }
      case ATAN2:
      {
        temp = std::make_shared<EBBinaryOpTermNode>(_left, _right, EBBinaryOpTermNode::NodeType::DIV);
        temp2 = std::make_shared<EBUnaryFuncTermNode>(temp, EBUnaryFuncTermNode::ATAN);
        derivative = temp2->takeDerivative();
        break;
      }
      case HYPOT:
      {
        temp = std::make_shared<EBBinaryOpTermNode>(_left, _right, EBBinaryOpTermNode::NodeType::DIV);
        temp2 = std::make_shared<EBBinaryOpTermNode>(temp, temp, EBBinaryOpTermNode::NodeType::MUL);
        temp3 = std::make_shared<EBBinaryOpTermNode>(std::make_shared<EBNumberNode<Real> >(1.0), temp2, EBBinaryOpTermNode::NodeType::ADD);
        temp4 = std::make_shared<EBUnaryFuncTermNode>(temp3, EBUnaryFuncTermNode::NodeType::SQRT);
        temp5 = std::make_shared<EBUnaryFuncTermNode>(_right, EBUnaryFuncTermNode::NodeType::ABS);
        temp6 = std::make_shared<EBBinaryOpTermNode>(temp4, temp5, EBBinaryOpTermNode::NodeType::MUL);
        derivative = temp6->takeDerivative();
        break;
      }
      case PLOG:
      {
        temp = std::make_shared<EBBinaryOpTermNode>(_left, _right, EBBinaryOpTermNode::NodeType::SUB);
        temp2 = std::make_shared<EBBinaryOpTermNode>(std::make_shared<EBNumberNode<Real> >(1.0), temp, EBBinaryOpTermNode::NodeType::ADD);
        temp3 = std::make_shared<EBBinaryOpTermNode>(std::make_shared<EBNumberNode<Real> >(1.0), temp, EBBinaryOpTermNode::NodeType::SUB);
        temp4 = std::make_shared<EBUnaryFuncTermNode>(temp2, EBUnaryFuncTermNode::NodeType::LOG);
        temp5 = std::make_shared<EBUnaryFuncTermNode>(temp3, EBUnaryFuncTermNode::NodeType::LOG);
        temp6 = std::make_shared<EBUnaryOpTermNode>(temp5, EBUnaryOpTermNode::NodeType::NEG);
        temp7 = std::make_shared<EBBinaryOpTermNode>(temp, std::make_shared<EBNumberNode<Real> >(0.0), EBBinaryOpTermNode::NodeType::LESS);
        temp8 = std::make_shared<EBTernaryFuncTermNode>(temp7, temp6, temp5, EBTernaryFuncTermNode::NodeType::CONDITIONAL);
        derivative = temp8->takeDerivative();
        break;
      }
    }
  }
  return derivative;
}

template <>
Real
ExpressionBuilder::EBMatPropNode<Real>::evaluate()
{
  return _mpd->value()[current_quad_point];
}

template <>
Real
ExpressionBuilder::EBMatPropNode<RealVectorValue>::evaluate()
{
  return _mpd->value()[current_quad_point](_comp1);
}

template <>
Real
ExpressionBuilder::EBMatPropNode<RankTwoTensor>::evaluate()
{
  return _mpd->value()[current_quad_point](_comp1, _comp2);
}

template <>
Real
ExpressionBuilder::EBCoupledVarNode<VariableValue>::evaluate()
{
  return (*_value)[current_quad_point];
}

template <>
Real
ExpressionBuilder::EBCoupledVarNode<VariableGradient>::evaluate()
{
  return ((*_value)[current_quad_point])(_comp1);
}

template <>
Real
ExpressionBuilder::EBCoupledVarNode<VariableSecond>::evaluate()
{
  return ((*_value)[current_quad_point])(_comp1, _comp2);
}

Real
ExpressionBuilder::EBUnaryFuncTermNode::evaluate()
{
  if(!isEvaluated)
  {
    switch (_type)
    {
      case SIN:
      {
        evaluated = sin(_subnode->evaluate());
        break;
      }
      case COS:
      {
        evaluated = cos(_subnode->evaluate());
        break;
      }
      case TAN:
      {
        evaluated = tan(_subnode->evaluate());
        break;
      }
      case ABS:
      {
        evaluated = abs(_subnode->evaluate());
        break;
      }
      case LOG:
      {
        evaluated = log(_subnode->evaluate());
        break;
      }
      case LOG2:
      {
        evaluated = log2(_subnode->evaluate());
        break;
      }
      case LOG10:
      {
        evaluated = log10(_subnode->evaluate());
        break;
      }
      case EXP:
      {
        evaluated = exp(_subnode->evaluate());
        break;
      }
      case COSH:
      {
        evaluated = cosh(_subnode->evaluate());
        break;
      }
      case SINH:
      {
        evaluated = sinh(_subnode->evaluate());
        break;
      }
      case SQRT:
      {
        evaluated = sqrt(_subnode->evaluate());
        break;
      }
      case ACOS:
      {
        evaluated = acos(_subnode->evaluate());
        break;
      }
      case ATAN:
      {
        evaluated = atan(_subnode->evaluate());
        break;
      }
    }
  }
  isEvaluated = true;
  return evaluated;
}

Real
ExpressionBuilder::EBBinaryOpTermNode::evaluate()
{
  if(!isEvaluated)
  {
    switch (_type)
    {
      case ADD:
      {
        evaluated = _left->evaluate() + _right->evaluate();
        break;
      }
      case SUB:
      {
        evaluated = _left->evaluate() - _right->evaluate();
        break;
      }
      case MUL:
      {
        evaluated = _left->evaluate() * _right->evaluate();
        break;
      }
      case DIV:
      {
        evaluated = _left->evaluate() / _right->evaluate();
        break;
      }
      case MOD:
      {
        evaluated = fmod(_left->evaluate(), _right->evaluate());
        break;
      }
      case POW:
      {
        evaluated = pow(_left->evaluate(), _right->evaluate());
        break;
      }
      case LESS:
      case GREATER:
      case LESSEQ:
      case GREATEREQ:
      case EQ:
      case NOTEQ:
      case AND:
      case OR:
        mooseError("Not a Real value");
    }
  }
  isEvaluated = true;
  return evaluated;
}

bool
ExpressionBuilder::EBBinaryOpTermNode::boolEvaluate()
{
  if(!isBoolEvaluated)
  {
    switch (_type)
    {
      case LESS:
      {
        boolEvaluated = _left->evaluate() < _right->evaluate();
        break;
      }
      case GREATER:
      {
        boolEvaluated = _left->evaluate() > _right->evaluate();
        break;
      }
      case LESSEQ:
      {
        boolEvaluated = _left->evaluate() <= _right->evaluate();
        break;
      }
      case GREATEREQ:
      {
        boolEvaluated = _left->evaluate() >= _right->evaluate();
        break;
      }
      case EQ:
      {
        boolEvaluated = _left->evaluate() == _right->evaluate();
        break;
      }
      case NOTEQ:
      {
        boolEvaluated = _left->evaluate() != _right->evaluate();
        break;
      }
      case AND:
      {
        boolEvaluated = _left->boolEvaluate() && _right->boolEvaluate();
        break;
      }
      case OR:
      {
        boolEvaluated = _left->boolEvaluate() || _right->boolEvaluate();
        break;
      }
      case SUB:
      case ADD:
      case MUL:
      case DIV:
      case MOD:
      case POW:
        mooseError("Not a boolean value");
    }
  }
  isBoolEvaluated = true;
  return boolEvaluated;
}

Real
ExpressionBuilder::EBBinaryFuncTermNode::evaluate()
{
  if(!isEvaluated)
  {
    Real left = _left->evaluate();
    Real right = _right->evaluate();
    switch (_type)
    {
      case MIN:
      {
        if(left < right)
          evaluated = left;
        else
          evaluated = right;
        break;
      }
      case MAX:
      {
        if(left > right)
          evaluated = left;
        else
          evaluated = right;
        break;
      }
      case ATAN2:
      {
        if(right > 0)
          evaluated = atan(left / right);
        else if(right < 0 && left >= 0)
          evaluated = atan(left / right) + pi;
        else if(right < 0 && left < 0)
          evaluated = atan(left / right) - pi;
        else if(right == 0 && left > 0)
          evaluated = pi / 2;
        else if(right == 0 && left < 0)
          evaluated = -pi / 2;
        else if(right == 0 && left == 0)
          evaluated = 0;
          //mooseError("atan2 is undefined at x = 0 and y = 0");
        break;
      }
      case HYPOT:
      {
        evaluated = abs(left) * sqrt(1.0 + right * right/(left * left));
        break;
      }
      case PLOG:
      {
        Real yprime = left - right;
        if(yprime >= 0)
          evaluated = log(1 + yprime);
        else
          evaluated = -log(1 - yprime);
        break;
      }
    }
  }
  isEvaluated = true;
  return evaluated;
}
