"""A module to handle the json output of MadNkLO runs
TODO DOC
"""
from math import sqrt
import json
from numbers import Number

class Result(object):
    """A class that contains a single MadNkLO results
    """
    def result_from_dict(self,result_dict):
        self.value = float(result_dict["Cross section"])
        self.uncertainty = float(result_dict["MC uncertainty"])
        self.name = str(result_dict["name"])
        self.timestamp = str(result_dict["timestamp"])
        self.unit = str(result_dict["unit"])
        self.order = str(result_dict["order"])

    def result_from_value_uncertainty(self,value,uncertainty):
        self.value = float(value)
        self.uncertainty = float(uncertainty)
        if uncertainty == 0.:
            self.name = "[{:1.1e}]".format(self.value)
        else:
            self.name = "[{:1.1e}+/-{:1.1e}]".format(self.value,self.uncertainty)
        self.timestamp = ""
        self.unit = "pb"
        self.order = ""

    def __init__(self,*input):
        """Get a dictionnary with the following keys:
        -TODO

        :param result_dict:
        :type result_dict:
        """
        #TODO Input sanitization
        if len(input) == 1 and isinstance(input[0], dict):
            self.result_from_dict(input[0])

        elif len(input) ==2 and all(isinstance(n, Number) for n in input):
            self.result_from_value_uncertainty(*input)

        elif len(input) ==1 and isinstance(input[0],Number):
            self.result_from_value_uncertainty(input[0],0)

        else:
            raise ValueError("Results expect either a dictionary from a JSON result file or two numbers (value,error) as argument")



    def __add__(self, other):
        """TODO DOC """
        if isinstance(other,Result):
            if self.unit != other.unit:
                print "WARNING YOU ARE ADDING APPLES AND BANANAS !!! ({}+{})".format(self.unit,other.unit)
            result = Result(
                {
                    "Cross section": self.value+other.value,
                    "MC uncertainty": sqrt(self.uncertainty**2+other.uncertainty**2),
                    "name": self.name+"+"+other.name,
                    "order": self.order+"+"+other.order,
                    "timestamp": self.timestamp+"+"+other.timestamp,
                    "unit": self.unit
                }
            )
            return result
        elif isinstance(other,Number):
            result = self+Result(other,0.)
            return result

    def __radd__(self, other):
        return self.__add__(other)


    def __sub__(self, other):
        """TODO DOC """
        if isinstance(other,Result):
            if self.unit != other.unit:
                print "WARNING YOU ARE SUBTRACTING APPLES AND BANANAS !!! ({}+{})".format(self.unit,other.unit)

            result = Result(
                {
                    "Cross section": self.value-other.value,
                    "MC uncertainty": sqrt(self.uncertainty**2+other.uncertainty**2),
                    "name": self.name+"-"+other.name,
                    "order": self.order+"-"+other.order,
                    "timestamp": self.timestamp+"-"+other.timestamp,
                    "unit": self.unit
                }
            )
            return result
        elif isinstance(other, Number):
            result = self - Result(other, 0.)
            return result

    def __rsub__(self, other):
        """TODO DOC"""
        if isinstance(other, Number):
            result = Result(other, 0.)-self
            return result

    def __mul__(self, other):
        """TODO DOC"""
        if isinstance(other,Result):
            if self.unit != other.unit:
                print "WARNING YOU ARE SUBTRACTING APPLES AND BANANAS !!! ({}+{})".format(self.unit,other.unit)

            E1 = self.value
            E2 = other.value
            S1 = self.uncertainty**2
            S2 = other.uncertainty**2
            # Assuming normal variables,
            # Formula: Var(X*Y) = Var(X)*Var(Y) + Var(X)*E(Y)**2 + Var(Y)*E(X)**2
            # There should be a check that the higher moments are negligble

            result = Result(
                {
                    "Cross section": E1*E2,
                    "MC uncertainty": sqrt(S1*S2+S1*E2**2+S2*E1**2),
                    "name": self.name+"*"+other.name,
                    "order": self.order+"*"+other.order,
                    "timestamp": self.timestamp+"*"+other.timestamp,
                    "unit": self.unit
                }
            )
            return result
        elif isinstance(other, Number):
            result = Result(other, 0.)*self
            return result

    def __rmul__(self, other):
        """TODO DOC"""
        return self*other

    def is_zero(self,n_sigma=2):
        return abs(self.value) < abs(n_sigma*self.uncertainty)

    def __str__(self):
        return "{:4.4e} +/- {:4.4e}".format(self.value,self.uncertainty)

    def __repr__(self):
        return ('{} {}='.format(type(self),self.name))+self.__str__()


class ResultDict(object):
    """A dictionnary of Results initialized from a JSON file"""
    def __init__(self,json_path=None):
        self._dict = {}
        if json_path!=None:
            #TODO add exception handling
            with open(json_path,"r") as my_json_file:
                json_dict = json.loads(my_json_file.read())
            for key in json_dict:
                self._dict[key]=Result(json_dict[key])

    def __getitem__(self, item):
        return self._dict[item]

    def __setitem__(self, key, value):
        self._dict[key] = value

    def __str__(self):
        stringDict = ""
        for k in self._dict.keys():
            stringDict+="{}: {}\n".format(k,str(self[k]))
        return stringDict

##########################
# Specific objects
##########################

zero_pb = Result(0)