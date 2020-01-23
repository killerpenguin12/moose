//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "gtest/gtest.h"

#include "ExpressionBuilder.h"

class ExpressionBuilderMatrixTest : public ::testing::Test, public ExpressionBuilder
{
};

TEST_F(ExpressionBuilderMatrixTest, test)
{
  EBTensor test({{1,2,3},{4,5,6},{7,8,9}});
  std::cout << test;
  std::cout << test + test;
  std::cout << test - test;
  std::cout << test*test;
  std::cout << test*2;
  std::cout << 2*test;
  test.transpose();
  std::cout << test;
}
