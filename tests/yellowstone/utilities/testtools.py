#==============================================================================
#
#  TestTools
#
#  This is a collection of functions that are useful for running the PyReshaper
#  tests on the Yellowstone compute system.
#
#==============================================================================


class TestDatabase(object):

    def __init__(self, filename=None):
        # See if there is a user-defined testinfo file, otherwise look for
        # default
        abspath = ''
        if filename:
            abspath = os.path.abspath(filename)
        else:
            this_dir = os.path.dirname(__file__)
            abspath = os.path.join(this_dir, 'testinfo.json')

        # Try opening and reading the testinfo file
        self._database = {}
        try:
            dbfile = open(filename, 'r')
            self._database = dict(json.load(dbfile))
            dbfile.close()
        except:
            err_msg = 'Problem reading and parsing test info file: ' \
                + str(abspath)
            raise ValueError(err_msg)

    def list_tests(self):
        print 'Tests found in the Test Database are:'
        print
        for test_name in self._database:
            print '   ' + str(test_name)
