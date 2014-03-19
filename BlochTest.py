import unittest 

from bloch_simulator import bloch

class BlochTest(unittest.TestCase):
    """
    Tests basic functionality of Python
    Bloch simulator. All test results are based 
    on results of Matlab function. Any errors in that
    evaluation will be replicated here. 
    """

    input_args = ( (),
                   (),
                   ()
                 )

    expected_output = ( (),
                        (),
                        ()
                      )
    
    def test_simple_bloch_values(self):
        """
        Series of input and output pairs from Matlab's calls to function.
        """
        for val in range(len(input_args)):
            self.assertEqual(expected_output[val], bloch(*input_args))

if __name__ == "__main__":
    unittest.main()
