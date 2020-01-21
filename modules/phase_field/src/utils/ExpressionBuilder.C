//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ExpressionBuilder.h"

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
  const char * name[] = {"sin", "cos", "tan", "abs", "log", "log2", "log10", "exp", "sinh", "cosh"};
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
  const char * name[] = {"+", "-", "*", "/", "%", "^", "<", ">", "<=", ">=", "=", "!="};
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
  if (narg != _eval_arguments.size())
    mooseError("EBFunction is used wth a different number of arguments than it was defined with.");

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
    return ExpressionBuilder::EBTerm(new ExpressionBuilder::EBUnaryFuncTermNode(                   \
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

#define BINARY_FUNC_IMPLEMENT(op, OP)                                                              \
  ExpressionBuilder::EBTerm op(const ExpressionBuilder::EBTerm & left,                             \
                               const ExpressionBuilder::EBTerm & right)                            \
  {                                                                                                \
    mooseAssert(left._root != NULL,                                                                \
                "Empty term provided as first argument of function " #op "()");                    \
    mooseAssert(right._root != NULL,                                                               \
                "Empty term provided as second argument of function " #op "()");                   \
    return ExpressionBuilder::EBTerm(new ExpressionBuilder::EBBinaryFuncTermNode(                  \
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
  return ExpressionBuilder::EBTerm(new ExpressionBuilder::EBBinaryOpTermNode(
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
    return ExpressionBuilder::EBTerm(new ExpressionBuilder::EBTernaryFuncTermNode(                 \
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
    EBTermNode * replace = rules[i]->apply(_subnode);
    if (replace != NULL)
    {
      delete _subnode;
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
    EBTermNode * replace = rules[i]->apply(_left);
    if (replace != NULL)
    {
      delete _left;
      _left = replace;
      success = 1;
      break;
    }
  }

  if (success == 0)
    success += _left->substitute(rules);

  for (unsigned int i = 0; i < nrule; ++i)
  {
    EBTermNode * replace = rules[i]->apply(_right);
    if (replace != NULL)
    {
      delete _right;
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
  EBTermNode * replace;

  for (unsigned int i = 0; i < nrule; ++i)
  {
    replace = rules[i]->apply(_left);
    if (replace)
    {
      delete _left;
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
      delete _middle;
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
      delete _right;
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
    EBTermNode * replace = rules[i]->apply(_root);
    if (replace != NULL)
    {
      delete _root;
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
  const EBSymbolNode * find_root = dynamic_cast<const EBSymbolNode *>(find.getRoot());
  if (find_root == NULL)
    mooseError("Function arguments must be pure symbols.");
  _find = find_root->stringify();

  // the term we want to substitute with
  if (replace.getRoot() != NULL)
    _replace = replace.cloneRoot();
  else
    mooseError("Trying to substitute in an empty term for ", _find);
}

ExpressionBuilder::EBTermNode *
ExpressionBuilder::EBTermSubstitution::substitute(const EBSymbolNode & node) const
{
  if (node.stringify() == _find)
    return _replace->clone();
  else
    return NULL;
}

ExpressionBuilder::EBTermNode *
ExpressionBuilder::EBLogPlogSubstitution::substitute(const EBUnaryFuncTermNode & node) const
{
  if (node._type == EBUnaryFuncTermNode::LOG)
    return new EBBinaryFuncTermNode(
        node.getSubnode()->clone(), _epsilon->clone(), EBBinaryFuncTermNode::PLOG);
  else
    return NULL;
}

template <typename U>
void
ExpressionBuilder::EBTensor::initialize(std::initializer_list<U> u, std::size_t _depth)
{
  // outer recursions
  checkDimensions(u.size(), _depth);
  for (auto t : u)
    initialize(t, _depth + 1);
}

void
ExpressionBuilder::EBTensor::initialize(std::initializer_list<int> t, std::size_t _depth)
{
  // innermost recursion
  checkDimensions(t.size(), _depth);
  for (const int & v : t)
    _data.push_back(EBTerm(v));
}

void
ExpressionBuilder::EBTensor::initialize(std::initializer_list<Real> t, std::size_t _depth)
{
  // innermost recursion
  checkDimensions(t.size(), _depth);
  for (const Real & v : t)
    _data.push_back(EBTerm(v));
}

void
ExpressionBuilder::EBTensor::initialize(std::initializer_list<std::string> t, std::size_t _depth)
{
  // innermost recursion
  checkDimensions(t.size(), _depth);
  for (const std::string & v : t)
    _data.push_back(EBTerm(v.c_str()));
}

void
ExpressionBuilder::EBTensor::checkDimensions(std::size_t size, std::size_t _depth)
{
  // check size
  if (_shape.size() <= _depth)
  {
    _access_data.resize(_depth + 1);
    unsigned int index_weight = 1;
    if (_access_data.size() == 1)
      index_weight = 0;
    else
      for (unsigned int i = 0; i < _shape.size(); ++i)
        index_weight *= _shape[i];
    _access_data[_depth] = index_weight;
    _shape.resize(_depth + 1);
    _shape[_depth] = size;
  }
  else
  {
    if (_shape[_depth] != size)
      throw std::runtime_error("Inconsistent initializer list size");
  }
}

void
ExpressionBuilder::EBTensor::printDebug()
{
  std::cout << "shape:";
  for (auto s : _shape)
    std::cout << ' ' << s;
  std::cout << '\n';

  std::cout << "data: ";
  for (auto d : _data)
    std::cout << ' ' << d;
  std::cout << '\n';
}

// operators

ExpressionBuilder::EBTerm
ExpressionBuilder::EBTensor::operator()(unsigned int index, unsigned int index_rest...)
{
  return (*this)(std::vector<unsigned int>(1, index), index_rest);
}

ExpressionBuilder::EBTerm
ExpressionBuilder::EBTensor::
operator()(std::vector<unsigned int> index_caught, unsigned int index, unsigned int index_rest...)
{
  index_caught.push_back(index);
  return (*this)(index_caught, index_rest);
}

ExpressionBuilder::EBTerm
ExpressionBuilder::EBTensor::operator()(std::vector<unsigned int> index_caught, unsigned int index)
{
  index_caught.push_back(index);
  unsigned int location = 0;
  for (unsigned int i = 0; i < _access_data.size(); ++i)
    location += _access_data[i] * index_caught[i];
  return _data[location];
}

ExpressionBuilder::EBTensor &
ExpressionBuilder::EBTensor::operator+(const EBTensor & rhs)
{
  if (this->_shape != rhs._shape)
    mooseError("Improper shape for addition.");
  EBTensor * result = new EBTensor();
  std::vector<EBTerm> new_data(this->_data.size());
  for (unsigned int i = 0; i < this->_data.size(); ++i)
    new_data[i] = this->_data[i] + rhs._data[i];
  result->_data = new_data;
  result->_shape = this->_shape;
  result->_access_data = this->_access_data;
  return *result;
}

ExpressionBuilder::EBTensor &
ExpressionBuilder::EBTensor::operator-(const EBTensor & rhs)
{
  if (this->_shape != rhs._shape)
    mooseError("Improper shape for addition.");
  EBTensor * result = new EBTensor();
  std::vector<EBTerm> new_data(this->_data.size());
  for (unsigned int i = 0; i < this->_data.size(); ++i)
    new_data[i] = this->_data[i] - rhs._data[i];
  result->_data = new_data;
  result->_shape = this->_shape;
  result->_access_data = this->_access_data;
  return *result;
}

ExpressionBuilder::EBTensor & ExpressionBuilder::EBTensor::operator*(const EBTensor & rhs)
{
  if (this->_shape.back() != rhs._shape[0])
    mooseError("Improper shape for multiplication.");
  EBTensor * result = new EBTensor();
  result->_shape = std::vector<long unsigned int>(this->_shape.begin(), this->_shape.end() - 1);
  result->_shape.insert(result->_shape.end(), rhs._shape.begin() + 1, rhs._shape.end());
  unsigned int data_size = 1;
  for (unsigned int i = 0; i < result->_shape[i]; ++i)
    data_size *= result->_shape[i];

  std::vector<EBTerm> new_data(data_size, EBTerm(0.0));
  unsigned int divider = 1;
  for (unsigned int i = 1; i < rhs._shape.size(); ++i)
    divider *= rhs._shape[i];
  for (unsigned int i = 0; i < data_size; ++i)
    for (unsigned int j = 0; j < rhs._shape[0]; ++j)
      new_data[i] +=
          this->_data[i / divider * rhs._shape[0] + j] * rhs._data[i % divider + divider * j];
  result->_data = new_data;
  return *result;
}

ExpressionBuilder::EBTensor &
ExpressionBuilder::EBTensor::contractMult(const EBTensor & lhs,
                                          const EBTensor & rhs,
                                          unsigned int contractor)
{
  for (unsigned int i = 0; i < contractor; ++i)
    if (this->_shape[this->_shape.size() - 1 - i] != rhs._shape[i])
      mooseError("Improper shape for multiplication.");
  EBTensor * result = new EBTensor();
  result->_shape =
      std::vector<long unsigned int>(this->_shape.begin(), this->_shape.end() - contractor);
  result->_shape.insert(result->_shape.end(), rhs._shape.begin() + contractor, rhs._shape.end());
  unsigned int data_size = 1;
  for (unsigned int i = 0; i < result->_shape[i]; ++i)
    data_size *= result->_shape[i];

  std::vector<EBTerm> new_data(data_size, EBTerm(0.0));
  unsigned int divider = 1;
  for (unsigned int i = contractor; i < rhs._shape.size(); ++i)
    divider *= rhs._shape[i];
  std::vector<long unsigned int> center(this->_shape.end() - contractor, this->_shape.end());
  EBTerm begin(0.0);
  for (unsigned int i = 0; i < data_size; ++i)
    new_data[i] = recurseMult(
        begin, i, divider, center, 0, std::vector<unsigned int>(center.size(), 0), lhs, rhs);
  result->_data = new_data;
  return *result;
}

ExpressionBuilder::EBTerm &
ExpressionBuilder::EBTensor::recurseMult(EBTerm & result_term,
                                         unsigned int index,
                                         unsigned int divider,
                                         std::vector<long unsigned int> center,
                                         unsigned int current,
                                         std::vector<unsigned int> current_vec,
                                         const EBTensor & lhs,
                                         const EBTensor & rhs)
{
  if (current < center.size())
    for (unsigned int i = 0; i < center[current]; ++i)
    {
      current_vec[current] = i;
      recurseMult(result_term, index, divider, center, current + 1, current_vec, lhs, rhs);
    }
  unsigned int lhs_index = 0;
  unsigned int multiplier = 1;
  unsigned int rhs_index = index % divider;
  for (unsigned int i = 1; i <= center.size(); ++i)
  {
    lhs_index += multiplier * current_vec[center.size() - i];
    rhs_index += multiplier * current_vec[center.size() - i] * divider;
    multiplier *= center[center.size() - 1];
  }
  lhs_index += multiplier * index / divider;
  EBTerm * new_result = new EBTerm(result_term + lhs._data[lhs_index] * rhs._data[rhs_index]);
  return *new_result;
}
