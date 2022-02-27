#Isaac Darling, May 2021
#revisit several calculus topics and implement them

#import arccos(x), sin(x), and sqrt(x) because their implementations are a series of operations involving inscrutible constants
from math import acos, sin, sqrt

#implements Euler's method
#WARNING: DO NOT USE PUBLICLY. this function is DANGEROUS due to the use of eval(), which can grant the end user remote code execution powers
def euler_method(point, eq, step, target):
    t, y = point
    if abs(target-t) <= 5e-5:#accounts for floating point rounding error
        return y
    return euler_method((t+step, y+eval(eq)*step), eq, step, target)

#warning: only defined in terms of real numbers
class Vector:
    #accepts a list of components and constructs a Vector() object with default value: <0, 0>
    def __init__(self, vals = None):
        vals = vals or [0, 0]
        self.components = list(vals)

    #effective class member for unit vector i
    @staticmethod
    def i():
        return Vector([1, 0, 0])

    #effective class member for unit vector j
    @staticmethod
    def j():
        return Vector([0, 1, 0])

    #effective class member for unit vector k
    @staticmethod
    def k():
        return Vector([0, 0, 1])

    #function for determining valid input for vector operations
    def match(self, other):
        if not isinstance(other, Vector):
            return False
        return len(self) == len(other)

    #returns magnitude of self
    def magnitude(self):
        return sqrt(sum([x*x for x in self.components]))

    #returns dot product of two vectors (including self)
    def dot(self, other):
        if not self.match(other):
            raise ArithmeticError("the dot product is only defined for vectors with an equal number of components")
        return sum([self.components[i]*other.components[i] for i in range(len(self))])

    #returns component of self in the direction of other
    def comp(self, other):
        if not isinstance(other, Vector):
            raise TypeError("two vectors are required to find the component of one in the directin of the other")
        mag = other.magnitude()
        if mag == 0:
            raise ValueError("the given vector is a zero vector")
        return self.dot(other)/mag

    #returns projection of self onto other
    def proj(self, other):
        if not isinstance(other, Vector):
            raise TypeError("two vectors are required to find the projection of one onto the other")
        return (self.dot(other)/other.dot(other)) * other

    #returns angle between self and other
    def angle(self, other):
        if not isinstance(other, Vector):
            raise TypeError("two vectors are required to find the angle between vectors")
        mag1 = self.magnitude()
        mag2 = other.magnitude()
        if mag1 == 0 or mag2 == 0:
            raise ValueError("at least one given vector is a zero vector")
        return acos(self.dot(other)/(mag1*mag2))

    #returns cross product of of two vectors (including self)
    def cross(self, other):
        if not (self.match(other) and len(self) == 3):
            raise ArithmeticError("the cross product is only well defined between three dimensional vectors")
        return ((self.components[1]*other.components[2]-self.components[2]*other.components[1])*Vector.i()
                - (self.components[0]*other.components[2]-self.components[2]*other.components[0])*Vector.j()
                + (self.components[0]*other.components[1]-self.components[1]*other.components[0])*Vector.k())

    #returns unit vector in the direction of self
    def unit(self):
        mag = self.magnitude()
        if mag == 0:
            raise ValueError("the given vector is a zero vector")
        return self * (1/mag)

    #returns normal vector of self and other
    def normal(self, other):
        if not isinstance(other, Vector):
            raise TypeError("two vectors are required to find their normal vector")
        cross_length = self.magnitude()*other.magnitude()*sin(self.angle(other))
        if cross_length == 0:
            raise ValueError("the given vectors are parallel, and therefore have no normal vector")
        return self.cross(other) * (1/cross_length)

    #construction for generating zero vectors of dimension n
    @staticmethod
    def zero(n):
        return Vector([0]*n)

    #returns volume of parallelepiped formed by u, v, and w
    @staticmethod
    def volume(u, v, w):
        if not (isinstance(u, Vector) and isinstance(v, Vector) and isinstance(w, Vector)):
            raise TypeError("this function requires three vectors as input")
        return abs(u.cross(v).dot(w))

    #returns vector addition between self and other
    def __add__(self, other):
        if not self.match(other):
            raise ArithmeticError("vectors can only be added to other vectors with an equal number of components")
        return Vector([self.components[i]+other.components[i] for i in range(len(self))])

    #returns vector subtraction between self and other
    def __sub__(self, other):
        if not self.match(other):
            raise ArithmeticError("vectors can only be subtracted from other vectors with an equal number of components")
        return Vector([self.components[i]-other.components[i] for i in range(len(self))])

    #returns scalar multiplication of self and other
    def __mul__(self, other):
        if not isinstance(other, (int, float)):
            raise TypeError("only defined for scalar multiplication using real numbers")
        return Vector([self.components[i]*other for i in range(len(self))])

    #calls __mul__ but allows for repositioning of the '*' operator
    def __rmul__(self, other):
        return self.__mul__(other)

    #returns the equality of self and other
    def __eq__(self, other):
        if not isinstance(other, Vector):
            return False
        return self.components == other.components

    #returns the number of components of self
    def __len__(self):
        return len(self.components)

    #returns string representation of self
    def __repr__(self):
        return f"<{str(self.components)[1:-1]}>"
