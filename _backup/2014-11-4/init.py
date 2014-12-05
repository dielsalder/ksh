class ClassOne:
    def __init__(self, a):
        self.a = a
        self.square()

    def square(self):
        self.a = self.a ** 2
        return self.a

a = ClassOne(2)
print a.a
