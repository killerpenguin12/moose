//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include <vector>
#include <ostream>
#include <sstream>
#include <iomanip>
#include <math.h>

#include "MooseError.h"
#include "libmesh/libmesh_common.h"
#include "RotationTensor.h"
#include "FunctionMaterialPropertyDescriptor.h"

using namespace libMesh;

/**
 * ExpressionBuilder adds an interface to derived classes that enables
 * convenient construction of FParser expressions through operator overloading.
 * It exposes the new types EBTerm and EBFunction
 * Variables used in your expressions are of type EBTerm. The following declares
 * three variables that can be used in an expression:
 *
 * EBTerm c1("c1"), c2("c3"), phi("phi");
 *
 * Declare a function 'G' and define it. Note the double bracket syntax '(())'':
 *
 * EBFunction G;
 * G((c1, c2, c3)) = c1 + 2 * c2 + 3 * pow(c3, 2);
 *
 * Performing a substitution is as easy as:
 * EBFunction H;
 * H((c1, c2)) = G((c1, c2, 1-c1-c2))
 *
 * Use the ```<<``` io operator to output functions or terms. Or use explicit or
 * implicit casts from EBFunction to std::string``` to pass a function to the
 * FParser Parse method. FParser variables are built using the ```args()``` method.
 *
 * FunctionParserADBase<Real> GParser;
 * GParser.Parse(G, G.args);
 */
class ExpressionBuilder
{
public:
  ExpressionBuilder(){};

  // forward delcarations
  class EBTerm;
  class EBTermNode;
  class EBFunction;
  class EBMatrix;
  class EBVector;
  class EBQuaternion;
  class EBSubstitutionRule;
  typedef std::vector<EBTerm> EBTermList;
  typedef std::vector<std::shared_ptr<EBTermNode> > EBTermNodeList;
  typedef std::vector<const EBSubstitutionRule *> EBSubstitutionRuleList;

  static std::shared_ptr<EBTermNode> derivative_of;
  static unsigned int current_quad_point;

  /// Base class for nodes in the expression tree
  class EBTermNode : public std::enable_shared_from_this<EBTermNode>
  {
  public:
    EBTermNode() : isSimplified(false), derivative(NULL), isEvaluated(false), isBoolEvaluated(false) {};
    virtual ~EBTermNode(){};

    virtual std::string stringify() const = 0;
    virtual unsigned int substitute(const EBSubstitutionRuleList & /*rule*/) { return 0; }
    virtual int precedence() const = 0;
    friend std::ostream & operator<<(std::ostream & os, const EBTermNode & node)
    {
      return os << node.stringify();
    }

    virtual std::shared_ptr<EBTermNode> simplify() = 0;
    virtual std::shared_ptr<EBTermNode> takeDerivative() = 0;
    virtual Real evaluate() = 0;
    virtual bool boolEvaluate() = 0;
    virtual void cleanUp()
    {
      simplified = NULL;
      isSimplified = false;
      derivative = NULL;
      isEvaluated = false;
      isBoolEvaluated = false;
    }

    std::shared_ptr<EBTermNode> simplified;
    bool isSimplified;
    std::shared_ptr<EBTermNode> derivative;
    bool isEvaluated;
    Real evaluated;
    bool isBoolEvaluated;
    bool boolEvaluated;
  };

  /// Template class for leaf nodes holding numbers in the expression tree
  template <typename T>
  class EBNumberNode : public EBTermNode
  {
    T _value;

  public:
    EBNumberNode(T value) : _value(value){};

    virtual std::string stringify() const;
    virtual int precedence() const { return 0; }

    virtual std::shared_ptr<EBTermNode> simplify()
    {
      if(_value == 0)
        return NULL;
      return this->shared_from_this();
    }
    virtual std::shared_ptr<EBTermNode> takeDerivative()
    {
      return std::make_shared<EBNumberNode<int> >(0);
    }

    virtual Real evaluate()
    {
      return double(_value);
    }
    virtual bool boolEvaluate()
    {
      mooseError("Not a defined type");
      return false;
    }
  };

  /// Template class for leaf nodes holding symbols (i.e. variables) in the expression tree
  class EBSymbolNode : public EBTermNode
  {
    std::string _symbol;

  public:
    EBSymbolNode(std::string symbol) : _symbol(symbol){};

    virtual std::string stringify() const;
    virtual int precedence() const { return 0; }
    virtual std::shared_ptr<EBTermNode> simplify()
    {
      return this->shared_from_this();
    }
    virtual std::shared_ptr<EBTermNode> takeDerivative()
    {
      mooseError("This is not defined for EBMaterials");
      return NULL;
    }
    virtual Real evaluate()
    {
      mooseError("Not a defined type");
      return 0;
    }
    virtual bool boolEvaluate()
    {
      mooseError("Not a defined type");
      return false;
    }

  };

  template <typename T2>
  class EBMatPropNode : public EBSymbolNode
  {
    FunctionMaterialPropertyDescriptor<T2> * _mpd;
    int _comp1;
    int _comp2;

  public:
    EBMatPropNode(const std::string & name, FunctionMaterialPropertyDescriptor<T2> * mpd, int comp1 = -1, int comp2 = -1)
      : EBSymbolNode(name), _mpd(mpd), _comp1(comp1), _comp2(comp2) {};

    virtual std::shared_ptr<EBTermNode> takeDerivative()
    {
      return std::make_shared<EBNumberNode<int> >(0);
    }
    virtual Real evaluate();
    virtual bool boolEvaluate()
    {
      mooseError("Not a defined type");
      return false;
    }
  };

  template <typename T2, bool = true>
  class EBCoupledVarNode : public EBSymbolNode
  {
    int _comp1;
    int _comp2;
    const T2 * _value;

  public:
    EBCoupledVarNode(const std::string & name, const T2 * value,
                     int comp1 = -1, int comp2 = -1)
      :  EBSymbolNode(name), _comp1(comp1), _comp2(comp2),
         _value(value){};

    virtual std::shared_ptr<EBTermNode> takeDerivative()
    {
      if(derivative_of == this->shared_from_this())
        return std::make_shared<EBNumberNode<Real> >(1.0);
      return std::make_shared<EBNumberNode<int> >(0);
    }
    virtual Real evaluate();
    virtual bool boolEvaluate()
    {
      mooseError("Not a defined type");
      return false;
    }
  };

  /**
   * Template class for leaf nodes holding anonymous IDs in the expression tree.
   * No such node must be left in the final expression that is serialized and passed to FParser
   */
  class EBTempIDNode : public EBTermNode
  {
    unsigned long _id;

  public:
    EBTempIDNode(unsigned int id) : _id(id){};

    virtual std::string stringify() const; // returns "[idnumber]"
    virtual int precedence() const { return 0; }

    virtual std::shared_ptr<EBTermNode> simplify()
    {
      return this->shared_from_this();
    }
    virtual std::shared_ptr<EBTermNode> takeDerivative()
    {
      mooseError("Anonymous nodes still present when ready to evaluate");
      return NULL;
    }
    virtual Real evaluate()
    {
      mooseError("Anonymous nodes still present when ready to evaluate");
      return 0;
    }
    virtual bool boolEvaluate()
    {
      mooseError("Not a defined type");
      return false;
    }
  };

  /// Base class for nodes with a single sub node (i.e. functions or operators taking one argument)
  class EBUnaryTermNode : public EBTermNode
  {
  public:
    EBUnaryTermNode(std::shared_ptr<EBTermNode> subnode) : _subnode(subnode){};
    virtual ~EBUnaryTermNode() {};

    virtual unsigned int substitute(const EBSubstitutionRuleList & rule);
    const std::shared_ptr<EBTermNode> getSubnode() const { return _subnode; }

    virtual void cleanUp()
    {
      simplified = NULL;
      isSimplified = false;
      derivative = NULL;
      isEvaluated = false;
      isBoolEvaluated = false;
      if(_subnode->isEvaluated || _subnode->isBoolEvaluated)
        _subnode->cleanUp();
    }

  protected:
    std::shared_ptr<EBTermNode> _subnode;
  };

  /// Node representing a function with two arguments
  class EBUnaryFuncTermNode : public EBUnaryTermNode
  {
  public:
    enum NodeType
    {
      SIN,
      COS,
      TAN,
      ABS,
      LOG,
      LOG2,
      LOG10,
      EXP,
      SINH,
      COSH,
      SQRT,
      ACOS,
      ATAN
    } _type;

    EBUnaryFuncTermNode(std::shared_ptr<EBTermNode> subnode, NodeType type)
      : EBUnaryTermNode(subnode), _type(type){};

    virtual std::string stringify() const;
    virtual int precedence() const { return 2; }

    virtual std::shared_ptr<EBTermNode> simplify();
    virtual std::shared_ptr<EBTermNode> takeDerivative();
    virtual Real evaluate();
    virtual bool boolEvaluate()
    {
      mooseError("Not a defined type");
      return false;
    }
  };

  /// Node representing a unary operator
  class EBUnaryOpTermNode : public EBUnaryTermNode
  {
  public:
    enum NodeType
    {
      NEG,
      LOGICNOT
    } _type;

    EBUnaryOpTermNode(std::shared_ptr<EBTermNode> subnode, NodeType type)
      : EBUnaryTermNode(subnode), _type(type){};

    virtual std::string stringify() const;
    virtual int precedence() const { return 3; }

    virtual std::shared_ptr<EBTermNode> simplify()
    {
      if(!isSimplified)
      {
        std::shared_ptr<EBTermNode> new_sub = _subnode->simplify();
        if(new_sub == NULL)
        {
          _subnode = std::make_shared<EBNumberNode<int> >(0);
          if(_type == NEG)
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

    virtual std::shared_ptr<EBTermNode> takeDerivative()
    {
      if(derivative == NULL)
      {
        switch (_type)
        {
          case NEG:
          {
            derivative = std::make_shared<EBUnaryOpTermNode>(_subnode->takeDerivative(), NodeType::NEG);
            break;
          }
          case LOGICNOT:
            mooseError("Logical Operator not in logical statement");
        }
      }
      return derivative;
    }

    virtual Real evaluate()
    {
      if(!isEvaluated)
      {
        switch (_type)
        {
          case NEG:
          {
            evaluated = -_subnode->evaluate();
            break;
          }
          case LOGICNOT:
            mooseError("not a Real value");
        }
      }
      isEvaluated = true;
      return evaluated;
    }
    virtual bool boolEvaluate()
    {
      if(!isBoolEvaluated)
      {
        switch (_type)
        {
          case LOGICNOT:
          {
            boolEvaluated = !_subnode->boolEvaluate();
            break;
          }
          case NEG:
            mooseError("not a boolean value");
        }
      }
      isBoolEvaluated = true;
      return boolEvaluated;
    }
  };

  /// Base class for nodes with two sub nodes (i.e. functions or operators taking two arguments)
  class EBBinaryTermNode : public EBTermNode
  {
  public:
    EBBinaryTermNode(std::shared_ptr<EBTermNode> left, std::shared_ptr<EBTermNode> right) : _left(left), _right(right){};
    virtual ~EBBinaryTermNode() {};

    virtual unsigned int substitute(const EBSubstitutionRuleList & rule);

    virtual void cleanUp()
    {
      simplified = NULL;
      isSimplified = false;
      derivative = NULL;
      isEvaluated = false;
      isBoolEvaluated = false;
      if(_left->isEvaluated || _left->isBoolEvaluated)
        _left->cleanUp();
      if(_right->isEvaluated || _right->isBoolEvaluated)
        _right->cleanUp();
    }

  protected:
    std::shared_ptr<EBTermNode> _left;
    std::shared_ptr<EBTermNode> _right;
  };

  /// Node representing a binary operator
  class EBBinaryOpTermNode : public EBBinaryTermNode
  {
  public:
    enum NodeType
    {
      ADD,
      SUB,
      MUL,
      DIV,
      MOD,
      POW,
      LESS,
      GREATER,
      LESSEQ,
      GREATEREQ,
      EQ,
      NOTEQ,
      AND,
      OR
    };

    EBBinaryOpTermNode(std::shared_ptr<EBTermNode> left, std::shared_ptr<EBTermNode> right, NodeType type)
      : EBBinaryTermNode(left, right), _type(type){};

    virtual std::string stringify() const;
    virtual int precedence() const;

    virtual std::shared_ptr<EBTermNode> simplify();
    virtual std::shared_ptr<EBTermNode> takeDerivative();
    virtual Real evaluate();
    virtual bool boolEvaluate();

  protected:
    NodeType _type;
  };

  /// Node representing a function with two arguments
  class EBBinaryFuncTermNode : public EBBinaryTermNode
  {
  public:
    enum NodeType
    {
      MIN,
      MAX,
      ATAN2,
      HYPOT,
      PLOG
    } _type;

    EBBinaryFuncTermNode(std::shared_ptr<EBTermNode> left, std::shared_ptr<EBTermNode> right, NodeType type)
      : EBBinaryTermNode(left, right), _type(type){};

    virtual std::string stringify() const;
    virtual int precedence() const { return 2; }

    virtual std::shared_ptr<EBTermNode> simplify()
    {
      if(!isSimplified)
      {
        std::shared_ptr<EBTermNode> new_left = _left->simplify();
        std::shared_ptr<EBTermNode> new_right = _right->simplify();
        if(new_left == NULL)
          _left = std::make_shared<EBNumberNode<int> >(0);
        else if(_left != new_left)
          _left = new_left;
        if(new_right == NULL)
          _right = std::make_shared<EBNumberNode<int> >(0);
        else if(_right != new_right)
          _right = new_right;
        simplified = this->shared_from_this();
        isSimplified = true;
      }
      return simplified;
    }
    virtual std::shared_ptr<EBTermNode> takeDerivative();
    virtual Real evaluate();
    virtual bool boolEvaluate()
    {
      mooseError("not a Real value");
      return false;
    }
  };

  /// Base class for nodes with two sub nodes (i.e. functions or operators taking two arguments)
  class EBTernaryTermNode : public EBBinaryTermNode
  {
  public:
    EBTernaryTermNode(std::shared_ptr<EBTermNode> left, std::shared_ptr<EBTermNode> middle, std::shared_ptr<EBTermNode> right)
      : EBBinaryTermNode(left, right), _middle(middle){};
    virtual ~EBTernaryTermNode() {};

    virtual unsigned int substitute(const EBSubstitutionRuleList & rule);

    virtual void cleanUp()
    {
      simplified = NULL;
      isSimplified = false;
      derivative = NULL;
      isEvaluated = false;
      isBoolEvaluated = false;
      if(_left->isEvaluated || _left->isBoolEvaluated)
        _left->cleanUp();
      if(_middle->isEvaluated || _middle->isBoolEvaluated)
        _middle->cleanUp();
      if(_right->isEvaluated || _right->isBoolEvaluated)
        _right->cleanUp();
    }

  protected:
    std::shared_ptr<EBTermNode> _middle;
  };

  /// Node representing a function with three arguments
  class EBTernaryFuncTermNode : public EBTernaryTermNode
  {
  public:
    enum NodeType
    {
      CONDITIONAL
    } _type;

    EBTernaryFuncTermNode(std::shared_ptr<EBTermNode> left, std::shared_ptr<EBTermNode> middle, std::shared_ptr<EBTermNode> right, NodeType type)
      : EBTernaryTermNode(left, middle, right), _type(type){};

    virtual std::string stringify() const;
    virtual int precedence() const { return 2; }

    virtual std::shared_ptr<EBTermNode> simplify()
    {
      if(!isSimplified)
      {
        std::shared_ptr<EBTermNode> new_left = _left->simplify();
        std::shared_ptr<EBTermNode> new_right = _right->simplify();
        std::shared_ptr<EBTermNode> new_middle = _middle->simplify();
        if(new_left == NULL)
          _left = std::make_shared<EBNumberNode<int> >(0);
        else if(_left != new_left)
          _left = new_left;
        if(new_middle == NULL)
          _middle = std::make_shared<EBNumberNode<int> >(0);
        else if(_middle != new_middle)
          _middle = new_middle;
        if(new_right == NULL)
          _right = std::make_shared<EBNumberNode<int> >(0);
        else if(_right != new_right)
          _right = new_right;
        simplified = this->shared_from_this();
        isSimplified = true;
      }
      return simplified;
    }
    virtual std::shared_ptr<EBTermNode> takeDerivative()
    {
      if(derivative == NULL)
      {
        switch (_type)
        {
          case CONDITIONAL:
            derivative = std::make_shared<EBTernaryFuncTermNode>(_left, _middle->takeDerivative(), _right->takeDerivative(), _type);
        }
      }
      return derivative;
    }

    virtual Real evaluate()
    {
      if(!isEvaluated)
      {
        if(_left->boolEvaluate())
          evaluated = _middle->evaluate();
        else
          evaluated = _right->evaluate();
      }
      return evaluated;
    }
    virtual bool boolEvaluate()
    {
      mooseError("Not a defined type");
      return false;
    }
  };

  /**
   * Substitution rule functor base class to perform flexible term substitutions
   */
  class EBSubstitutionRule
  {
  public:
    virtual std::shared_ptr<EBTermNode> apply(const std::shared_ptr<EBTermNode>) const = 0;
    virtual ~EBSubstitutionRule() {}
  };

  /**
   * Substitution rule base class that applies to nodes of type Node_T
   */
  template <class Node_T>
  class EBSubstitutionRuleTyped : public EBSubstitutionRule
  {
  public:
    virtual std::shared_ptr<EBTermNode> apply(const std::shared_ptr<EBTermNode>) const;

  protected:
    // on successful substitution this returns a new node to replace the old one, otherwise it
    // returns NULL
    virtual std::shared_ptr<EBTermNode> substitute(const Node_T &) const = 0;
  };

  /**
   * Generic Substitution rule to replace all occurences of a given symbol node
   * term with a user defined term. This is used by EBFunction.
   */
  class EBTermSubstitution : public EBSubstitutionRuleTyped<EBSymbolNode>
  {
  public:
    EBTermSubstitution(const EBTerm & find, const EBTerm & replace);
    virtual ~EBTermSubstitution() {};

  protected:
    virtual std::shared_ptr<EBTermNode> substitute(const EBSymbolNode &) const;
    std::string _find;
    std::shared_ptr<EBTermNode> _replace;
  };

  /**
   * Substitution rule to replace all occurences of log(x) with plog(x, epsilon)
   * with a user defined term for epsilon.
   */
  class EBLogPlogSubstitution : public EBSubstitutionRuleTyped<EBUnaryFuncTermNode>
  {
  public:
    EBLogPlogSubstitution(const EBTerm & epsilon) : _epsilon(epsilon.cloneRoot())
    {
      mooseAssert(_epsilon != NULL, "Epsilon must not be an empty term in EBLogPlogSubstitution");
    }
    virtual ~EBLogPlogSubstitution() {};

  protected:
    virtual std::shared_ptr<EBTermNode> substitute(const EBUnaryFuncTermNode &) const;
    std::shared_ptr<EBTermNode> _epsilon;
  };

  /**
   * User facing host object for an expression tree. Each EBTerm contains a _root
   * node pointer to an EBTermNode object. The _root pointer should never be NULL,
   * but it should be safe if it ever is. The default constructor assigns a
   * EBTempIDNode to _root with a unique ID.
   */
  class EBTerm
  {
  public:
    // the default constructor assigns a temporary id node to root we use the address of the
    // current EBTerm object as the ID. This could be problematic if we create and destroy terms,
    // but then we should not expect the substitution to do sane things anyways.
    EBTerm() : _root(std::make_shared<EBTempIDNode>(reinterpret_cast<unsigned long>(this))){};

    EBTerm(const EBTerm & term) : _root(term.cloneRoot()){};
    ~EBTerm() {};

    // construct a term from a node
    EBTerm(std::shared_ptr<EBTermNode> root) : _root(root){};

  public:
    // construct from number or string
    EBTerm(int number) : _root(std::make_shared<EBNumberNode<int> >(number)) {}
    EBTerm(Real number) : _root(std::make_shared<EBNumberNode<Real> >(number)) {}
    EBTerm(const char * symbol) : _root(std::make_shared<EBSymbolNode>(symbol)) {}

    typedef std::vector<EBTerm> EBTermVector;
    static EBTermVector CreateEBTermVector(const std::string & var_name, unsigned int _op_num);

    // concatenate terms to form a parameter list with (()) syntax (those need to be out-of-class!)
    friend EBTermList operator,(const ExpressionBuilder::EBTerm & larg,
                                const ExpressionBuilder::EBTerm & rarg);
    friend EBTermList operator,(const ExpressionBuilder::EBTerm & larg,
                                const ExpressionBuilder::EBTermList & rargs);
    friend EBTermList operator,(const ExpressionBuilder::EBTermList & largs,
                                const ExpressionBuilder::EBTerm & rarg);

    // dump term as FParser expression
    friend std::ostream & operator<<(std::ostream & os, const EBTerm & term);
    // cast into a string
    operator std::string() const { return _root->stringify(); }

    // assign a term
    EBTerm & operator=(const EBTerm & term)
    {
      _root = term.cloneRoot();
      return *this;
    }

    // perform a substitution (returns substituton count)
    unsigned int substitute(const EBSubstitutionRule & rule);
    unsigned int substitute(const EBSubstitutionRuleList & rules);

    const std::shared_ptr<EBTermNode> getRoot() const { return _root; }
    std::shared_ptr<EBTermNode> cloneRoot() const { return _root; }

    bool simplify()
    {
      _root = _root->simplify();
      if(_root == NULL)
      {
        _root = std::make_shared<EBNumberNode<int> >(0);
        return false;
      }
      return true;
    }

    EBTerm takeDerivative()
    {
      return _root->takeDerivative();
    }

    Real evaluate()
    {
      return _root->evaluate();
    }

    void cleanUp()
    {
      _root->cleanUp();
    }

  protected:
    std::shared_ptr<EBTermNode> _root;

  public:
/**
 * Unary operators
 */
#define UNARY_OP_IMPLEMENT(op, OP)                                                                 \
  EBTerm operator op() const                                                                       \
  {                                                                                                \
    mooseAssert(_root != NULL, "Empty term provided for unary operator " #op);                     \
    return EBTerm(std::make_shared<EBUnaryOpTermNode>(cloneRoot(), EBUnaryOpTermNode::OP));                      \
  }
    UNARY_OP_IMPLEMENT(-, NEG)
    UNARY_OP_IMPLEMENT(!, LOGICNOT)

    /**
     * Unary functions
     */
    friend EBTerm sin(const EBTerm &);
    friend EBTerm cos(const EBTerm &);
    friend EBTerm tan(const EBTerm &);
    friend EBTerm abs(const EBTerm &);
    friend EBTerm log(const EBTerm &);
    friend EBTerm log2(const EBTerm &);
    friend EBTerm log10(const EBTerm &);
    friend EBTerm exp(const EBTerm &);
    friend EBTerm sinh(const EBTerm &);
    friend EBTerm cosh(const EBTerm &);
    friend EBTerm sqrt(const EBTerm &);
    friend EBTerm acos(const EBTerm &);
    friend EBTerm atan(const EBTerm &);

/*
 * Binary operators (including number,term operations)
 */
#define BINARY_OP_IMPLEMENT(op, OP)                                                                \
  EBTerm operator op(const EBTerm & term) const                                                    \
  {                                                                                                \
    mooseAssert(_root != NULL, "Empty term provided on left side of operator " #op);               \
    mooseAssert(term._root != NULL, "Empty term provided on right side of operator " #op);         \
    return EBTerm(std::make_shared<EBBinaryOpTermNode>(cloneRoot(), term.cloneRoot(), EBBinaryOpTermNode::OP));  \
  }                                                                                                \
  friend EBTerm operator op(int left, const EBTerm & right)                                        \
  {                                                                                                \
    mooseAssert(right._root != NULL, "Empty term provided on right side of operator " #op);        \
    return EBTerm(std::make_shared<EBBinaryOpTermNode>(                                            \
        std::make_shared<EBNumberNode<int> >(left), right.cloneRoot(), EBBinaryOpTermNode::OP));   \
  }                                                                                                \
  friend EBTerm operator op(Real left, const EBTerm & right)                                       \
  {                                                                                                \
    mooseAssert(right._root != NULL, "Empty term provided on right side of operator " #op);        \
    return EBTerm(std::make_shared<EBBinaryOpTermNode>(                                            \
        std::make_shared<EBNumberNode<Real> >(left), right.cloneRoot(), EBBinaryOpTermNode::OP));  \
  }                                                                                                \
  friend EBTerm operator op(const EBFunction & left, const EBTerm & right)                         \
  {                                                                                                \
    mooseAssert(EBTerm(left)._root != NULL, "Empty term provided on left side of operator " #op);  \
    mooseAssert(right._root != NULL, "Empty term provided on right side of operator " #op);        \
    return EBTerm(std::make_shared<EBBinaryOpTermNode>(                                            \
        EBTerm(left).cloneRoot(), right.cloneRoot(), EBBinaryOpTermNode::OP));                     \
  }                                                                                                \
  friend EBTerm operator op(const EBFunction & left, const EBFunction & right);                    \
  friend EBTerm operator op(int left, const EBFunction & right);                                   \
  friend EBTerm operator op(Real left, const EBFunction & right);
    BINARY_OP_IMPLEMENT(+, ADD)
    BINARY_OP_IMPLEMENT(-, SUB)
    BINARY_OP_IMPLEMENT(*, MUL)
    BINARY_OP_IMPLEMENT(/, DIV)
    BINARY_OP_IMPLEMENT(%, MOD)
    BINARY_OP_IMPLEMENT(<, LESS)
    BINARY_OP_IMPLEMENT(>, GREATER)
    BINARY_OP_IMPLEMENT(<=, LESSEQ)
    BINARY_OP_IMPLEMENT(>=, GREATEREQ)
    BINARY_OP_IMPLEMENT(==, EQ)
    BINARY_OP_IMPLEMENT(!=, NOTEQ)
    BINARY_OP_IMPLEMENT(&&, AND)
    BINARY_OP_IMPLEMENT(||,OR)

/*
 * Compound assignment operators
 */
#define BINARYCOMP_OP_IMPLEMENT(op, OP)                                                            \
  EBTerm & operator op(const EBTerm & term)                                                        \
  {                                                                                                \
    mooseAssert(_root != NULL, "Empty term provided on left side of operator " #op);               \
    mooseAssert(term._root != NULL, "Empty term provided on right side of operator " #op);         \
    if (std::dynamic_pointer_cast<EBTempIDNode>(_root))                                                       \
      mooseError("Using compound assignment operator on anonymous term. Set it to 0 first!");      \
    _root = std::make_shared<EBBinaryOpTermNode>(_root, term.cloneRoot(), EBBinaryOpTermNode::OP);               \
    return *this;                                                                                  \
  }
    BINARYCOMP_OP_IMPLEMENT(+=, ADD)
    BINARYCOMP_OP_IMPLEMENT(-=, SUB)
    BINARYCOMP_OP_IMPLEMENT(*=, MUL)
    BINARYCOMP_OP_IMPLEMENT(/=, DIV)
    BINARYCOMP_OP_IMPLEMENT(%=, MOD)

    /**
    * @{
     * Binary functions
     */
    friend EBTerm min(const EBTerm &, const EBTerm &);
    friend EBTerm max(const EBTerm &, const EBTerm &);
    friend EBTerm pow(const EBTerm &, const EBTerm &);
    template <typename T>
    friend EBTerm pow(const EBTerm &, T exponent);
    friend EBTerm atan2(const EBTerm &, const EBTerm &);
    friend EBTerm hypot(const EBTerm &, const EBTerm &);
    friend EBTerm plog(const EBTerm &, const EBTerm &);
    ///@}

    /**
     * Ternary functions
     */

    friend EBTerm conditional(const EBTerm &, const EBTerm &, const EBTerm &);
    /*
    friend EBVector conditional(const EBTerm &, const EBVector &, const EBVector &);
    friend EBMatrix conditional(const EBTerm &, const EBMatrix &, const EBMatrix &);
    friend EBQuaternion conditional(const EBTerm &, const EBQuaternion &, const EBQuaternion &);
    */
  };

  /// User facing host object for a function. This combines a term with an argument list.
  class EBFunction
  {
  public:
    EBFunction(){};

    /// @{
    /// set the temporary argument list which is either used for evaluation
    /// or committed to the argument list upon function definition (assignment)
    EBFunction & operator()(const EBTerm & arg);
    EBFunction & operator()(const EBTermList & args);
    /// @}

    /// @{
    /// convenience operators to allow single bracket syntax
    EBFunction & operator()(const EBTerm & a1, const EBTerm & a2) { return (*this)((a1, a2)); }
    EBFunction & operator()(const EBTerm & a1, const EBTerm & a2, const EBTerm & a3)
    {
      return (*this)((a1, a2, a3));
    }
    EBFunction &
    operator()(const EBTerm & a1, const EBTerm & a2, const EBTerm & a3, const EBTerm & a4)
    {
      return (*this)((a1, a2, a3, a4));
    }
    EBFunction & operator()(const EBTerm & a1,
                            const EBTerm & a2,
                            const EBTerm & a3,
                            const EBTerm & a4,
                            const EBTerm & a5)
    {
      return (*this)((a1, a2, a3, a4, a5));
    }
    EBFunction & operator()(const EBTerm & a1,
                            const EBTerm & a2,
                            const EBTerm & a3,
                            const EBTerm & a4,
                            const EBTerm & a5,
                            const EBTerm & a6)
    {
      return (*this)((a1, a2, a3, a4, a5, a6));
    }
    EBFunction & operator()(const EBTerm & a1,
                            const EBTerm & a2,
                            const EBTerm & a3,
                            const EBTerm & a4,
                            const EBTerm & a5,
                            const EBTerm & a6,
                            const EBTerm & a7)
    {
      return (*this)((a1, a2, a3, a4, a5, a6, a7));
    }
    EBFunction & operator()(const EBTerm & a1,
                            const EBTerm & a2,
                            const EBTerm & a3,
                            const EBTerm & a4,
                            const EBTerm & a5,
                            const EBTerm & a6,
                            const EBTerm & a7,
                            const EBTerm & a8)
    {
      return (*this)((a1, a2, a3, a4, a5, a6, a7, a8));
    }
    EBFunction & operator()(const EBTerm & a1,
                            const EBTerm & a2,
                            const EBTerm & a3,
                            const EBTerm & a4,
                            const EBTerm & a5,
                            const EBTerm & a6,
                            const EBTerm & a7,
                            const EBTerm & a8,
                            const EBTerm & a9)
    {
      return (*this)((a1, a2, a3, a4, a5, a6, a7, a8, a9));
    }
    /// @}

    /// cast an EBFunction into an EBTerm
    operator EBTerm() const;

    /// cast into a string (via the cast into a term above)
    operator std::string() const;

    /// @{
    /// function definition (assignment)
    EBFunction & operator=(const EBTerm &);
    EBFunction & operator=(const EBFunction &);
    /// @}

    /// get the list of arguments and check if they are all symbols
    std::string args();

    /// @{
    /// Unary operators on functions
    EBTerm operator-() { return -EBTerm(*this); }
    EBTerm operator!() { return !EBTerm(*this); }
    /// @}

    // perform a substitution (returns substituton count)
    unsigned int substitute(const EBSubstitutionRule & rule);
    unsigned int substitute(const EBSubstitutionRuleList & rules);

  protected:
    /// argument list the function is declared with
    EBTermList _arguments;
    /// argument list passed in when evaluating the function
    EBTermList _eval_arguments;

    // underlying term that the _eval_arguments are substituted in
    EBTerm _term;
  };

  class EBMatrix
  {
  public:
    EBMatrix();
    EBMatrix(std::vector<std::vector<EBTerm> > FunctionMatrix);
    EBMatrix(unsigned int i, unsigned int j);
    EBMatrix(const RealTensorValue & rhs);

    operator EBQuaternion() const;
    //operator MaterialPropertyVariable<RankTwoTensor>() const;
    operator std::vector<std::vector<std::string> >() const;
    friend EBMatrix & operator*(const EBMatrix & lhs, const EBMatrix & rhs);
    friend EBMatrix & operator*(const Real & lhs, const EBMatrix & rhs);
    friend EBMatrix & operator*(const EBMatrix & lhs, const Real & rhs);
    friend EBMatrix & operator*(const EBTerm & lhs, const EBMatrix & rhs);
    friend EBMatrix & operator*(const EBMatrix & lhs, const EBTerm & rhs);
    friend EBMatrix & operator/(const EBMatrix & lhs, const EBTerm & rhs);
    friend EBMatrix & operator/(const EBMatrix & lhs, const Real & rhs);
    friend EBMatrix & operator+(const EBMatrix & lhs, const EBMatrix & rhs);
    friend EBMatrix & operator-(const EBMatrix & lhs, const EBMatrix & rhs);
    EBMatrix & operator+=(const EBMatrix & rhs);
    EBMatrix & operator-();
    std::vector<EBTerm> & operator[](unsigned int i);
    const std::vector<EBTerm> & operator[](unsigned int i) const;
    EBMatrix transpose();

    unsigned int rowNum() const;
    unsigned int colNum() const;
    void setSize(const unsigned int i, const unsigned int j); //Erases the entire matrix
    static EBMatrix rotVec1ToVec2(EBVector & Vec1, EBVector & Vec2);
    static EBMatrix rotVecToZ(EBVector & Vec);

    void checkMultSize(const EBMatrix & lhs, const EBMatrix & rhs);
    static void checkAddSize(const EBMatrix & lhs, const EBMatrix & rhs);

    void simplify()
    {
      for(unsigned int i = 0; i < _rowNum; ++i)
        for(unsigned int j = 0; j < _colNum; ++j)
          FunctionMatrix[i][j].simplify();
    }

  private:
    void checkSize();

    std::vector<std::vector<EBTerm> > FunctionMatrix; // Row then column
    unsigned int _rowNum;
    unsigned int _colNum;
  };

  class EBVector //3D vectors only, else use EBMatrix
  {
  public:
    EBVector();
    EBVector(Real i, Real j, Real k)
    {
      FunctionVector.push_back(EBTerm(i));
      FunctionVector.push_back(EBTerm(j));
      FunctionVector.push_back(EBTerm(k));
    }
    EBVector(EBTerm i, EBTerm j, EBTerm k)
    {
      FunctionVector.push_back(i);
      FunctionVector.push_back(j);
      FunctionVector.push_back(k);
    }
    EBVector(std::vector<EBTerm> FunctionVector);
    EBVector(std::vector<std::string> FunctionNameVector);

    typedef std::vector<EBVector> EBVectorVector;
    static EBVectorVector CreateEBVectorVector(const std::string & var_name, unsigned int _op_num);

    EBVector & operator=(const RealVectorValue & rhs);
    operator std::vector<std::string>() const;
    friend EBTerm & operator*(const EBVector & lhs, const EBVector & rhs); // This defines a dot product
    friend EBVector & operator*(const EBVector & lhs, const Real & rhs);
    friend EBVector & operator*(const Real & lhs, const EBVector & rhs);
    friend EBVector & operator*(const EBTerm & lhs, const EBVector & rhs);
    friend EBVector & operator*(const EBVector & lhs, const EBTerm & rhs);
    friend EBVector & operator*(const EBVector & lhs, const EBMatrix & rhs); // We assume the vector is 1 x 3 here
    friend EBVector & operator*(const EBMatrix & lhs, const EBVector & rhs); // We assume the vector is 3 x 1 here
    friend EBVector & operator/(const EBVector & lhs, const Real & rhs);
    friend EBVector & operator/(const EBVector & lhs, const EBTerm & rhs);
    friend EBVector & operator+(const EBVector & lhs, const EBVector & rhs);
    friend EBVector & operator-(const EBVector & lhs, const EBVector & rhs);
    EBVector & operator+=(const EBVector & rhs);
    EBVector & operator-();
    EBTerm & operator[](unsigned int i);
    const EBTerm & operator[](unsigned int i) const;
    //operator EBMatrix();

    static EBVector cross(const EBVector & lhs, const EBVector & rhs);
    EBTerm norm();
    void push_back(EBTerm term);

    void simplify()
    {
      for(unsigned int i = 0; i < 3; ++i)
        FunctionVector[i].simplify();
    }

  private:
    void checkSize(std::vector<std::string> FunctionVector);
    void checkSize(std::vector<EBTerm> FunctionVector);

    std::vector<EBTerm> FunctionVector;
  };

  class EBQuaternion
  {
  public:
    EBQuaternion() {
      FunctionQuat = std::vector<EBTerm>(4);
    };
    EBQuaternion(std::vector<EBTerm> FunctionQuat);

    friend EBQuaternion & operator*(const EBQuaternion & lhs, const EBQuaternion & rhs);
    friend EBQuaternion & operator*(const EBQuaternion & lhs, const Real & rhs);
    friend EBQuaternion & operator*(const Real & lhs, const EBQuaternion & rhs);
    friend EBQuaternion & operator*(const EBQuaternion & lhs, const EBTerm & rhs);
    friend EBQuaternion & operator*(const EBTerm & lhs, const EBQuaternion & rhs);
    friend EBQuaternion & operator/(const EBQuaternion & lhs, const Real & rhs);
    friend EBQuaternion & operator/(const EBQuaternion & lhs, const EBTerm & rhs);
    friend EBQuaternion & operator+(const EBQuaternion & lhs, const EBQuaternion & rhs);
    friend EBQuaternion & operator-(const EBQuaternion & lhs, const EBQuaternion & rhs);
    EBQuaternion & operator-(const EBQuaternion & rhs);
    EBTerm & operator[](unsigned int i);
    const EBTerm & operator[](unsigned int i) const;
    operator EBMatrix() const;

    EBTerm norm();

  private:
    void checkSize(std::vector<EBTerm> FunctionQuat);

    std::vector<EBTerm> FunctionQuat;
  };

/*
 * Binary operators
 */
#define BINARYFUNC_OP_IMPLEMENT(op, OP)                                                            \
  friend EBTerm operator op(const EBFunction & left, const EBFunction & right)                     \
  {                                                                                                \
    mooseAssert(EBTerm(left)._root != NULL, "Empty term provided on left side of operator " #op);  \
    mooseAssert(EBTerm(right)._root != NULL,                                                       \
                "Empty term provided on right side of operator " #op);                             \
    return EBTerm(std::make_shared<EBBinaryOpTermNode>(                                                          \
        EBTerm(left).cloneRoot(), EBTerm(right).cloneRoot(), EBBinaryOpTermNode::OP));             \
  }                                                                                                \
  friend EBTerm operator op(int left, const EBFunction & right)                                    \
  {                                                                                                \
    mooseAssert(EBTerm(right)._root != NULL,                                                       \
                "Empty term provided on right side of operator " #op);                             \
    return EBTerm(std::make_shared<EBBinaryOpTermNode>(                                                          \
        std::make_shared<EBNumberNode<int> >(left), EBTerm(right).cloneRoot(), EBBinaryOpTermNode::OP));          \
  }                                                                                                \
  friend EBTerm operator op(Real left, const EBFunction & right)                                   \
  {                                                                                                \
    mooseAssert(EBTerm(right)._root != NULL,                                                       \
                "Empty term provided on right side of operator " #op);                             \
    return EBTerm(std::make_shared<EBBinaryOpTermNode>(                                                          \
        std::make_shared<EBNumberNode<Real> >(left), EBTerm(right).cloneRoot(), EBBinaryOpTermNode::OP));         \
  }
  BINARYFUNC_OP_IMPLEMENT(+, ADD)
  BINARYFUNC_OP_IMPLEMENT(-, SUB)
  BINARYFUNC_OP_IMPLEMENT(*, MUL)
  BINARYFUNC_OP_IMPLEMENT(/, DIV)
  BINARYFUNC_OP_IMPLEMENT(%, MOD)
  BINARYFUNC_OP_IMPLEMENT(<, LESS)
  BINARYFUNC_OP_IMPLEMENT(>, GREATER)
  BINARYFUNC_OP_IMPLEMENT(<=, LESSEQ)
  BINARYFUNC_OP_IMPLEMENT(>=, GREATEREQ)
  BINARYFUNC_OP_IMPLEMENT(==, EQ)
  BINARYFUNC_OP_IMPLEMENT(!=, NOTEQ)
  BINARYFUNC_OP_IMPLEMENT(&&, AND)
  BINARYFUNC_OP_IMPLEMENT(||, OR)
};

// convenience function for numeric exponent
template <typename T>
ExpressionBuilder::EBTerm
pow(const ExpressionBuilder::EBTerm & left, T exponent)
{
  return ExpressionBuilder::EBTerm(
      std::make_shared<ExpressionBuilder::EBBinaryOpTermNode>(left.cloneRoot(),
                                                std::make_shared<ExpressionBuilder::EBNumberNode<T> >(exponent),
                                                ExpressionBuilder::EBBinaryOpTermNode::POW));
}

// convert a number node into a string
template <typename T>
std::string
ExpressionBuilder::EBNumberNode<T>::stringify() const
{
  std::ostringstream s;
  s << std::setprecision(12) << _value;
  return s.str();
}

template <class Node_T>
std::shared_ptr<ExpressionBuilder::EBTermNode>
ExpressionBuilder::EBSubstitutionRuleTyped<Node_T>::apply(
    const std::shared_ptr<ExpressionBuilder::EBTermNode> node) const
{
  const std::shared_ptr<Node_T> match_node = std::dynamic_pointer_cast<Node_T>(node);
  if (match_node == NULL)
    return NULL;
  else
    return substitute(*match_node);
}
