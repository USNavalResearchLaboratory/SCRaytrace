/*! \file boostTestTest.cpp 
 * \brief Test of the Boost unit test framework
 *
 *  This is just to make sure that it works...
 */

#define BOOST_TEST_MODULE example
// #include <boost/test/included/unit_test.hpp>

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE( free_test_function )
/* Compare with void free_test_function() */
{
  BOOST_TEST( true /* test assertion */ );
}
