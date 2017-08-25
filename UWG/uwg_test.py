class UWG_Test(object):
    def __init__(self,test_name,run_test=True):
        self.fail = 0
        self.success = 0
        self.total = self._get_total()
        self.test_history = ""
        self.test_name = test_name
        self.run_test = run_test
    def __repr__(self):
        return "{a} successful and {b} failed tests".format(a=self.success,b=self.fail)
    def test_results(self):
        if not self.run_test:
            return "\n-----TEST '" + self.test_name + "' NOT RUN-----"
        else:
            return "\n-----TEST '" + self.test_name + "' RESULT-----\n" + self.test_history + self.__repr__() +\
             "\n-----TEST '" + self.test_name + "' RESULTS-----"
    def _get_total(self):
        return self.fail + self.success
    def test_equality(self,a,b,toggle=True):
        if not self.run_test: return None
        if a == b:
            s = "test_equality: {y} == {z} success\n".format(y=a,z=b)
            self.success += 1
        else:
            s = "test_equality: {y} != {z} fail\n".format(y=a,z=b)
            self.fail += 1
        if toggle:
            self.test_history += s
    def test_equality_tol(self,a,b,toggle=True):
        if not self.run_test: return None
        tol = 0.001
        if abs(a-b) < 0.001:
            s = "test_equality_tol: {y} == {z} success\n".format(y=a,z=b)
            self.success += 1
        else:
            s = "test_equality_tol: {y} != {z} fail\n".format(y=a,z=b)
            self.fail += 1
        if toggle:
            self.test_history += s
    def test_in_string(self,a,b,toggle=True):
        if not self.run_test: return None
        if type(a)!=type("") or type(b)!=type(""):
            s = "test_in_string: {y} or {z} not a string\n".format(y=b,z=a)
            self.fail += 1
        else:
            if b in a:
                s = "test_in_string:: {y} in {z} success\n".format(y=b,z=a)
                self.success += 1
            else:
                s = "test_in_string: {y} in {z} fail\n".format(y=b,z=a)
                self.fail += 1
        if toggle:
            self.test_history += s
