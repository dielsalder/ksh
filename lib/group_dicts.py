import operator
import itertools

lst = [{'A':12,'B':32,'ID':333},{'Z':32,'C':43,'ID':111},{'D':43,'J':31,'ID':222},{'a':32,'b':31,'ID':222},{'D':43,'ID':333},{'a':89,'d':31,'ID':222},{'C':83,'ID':111}]

keyfunc = operator.itemgetter("ID")
desired_list = [list(grp) for key, grp in itertools.groupby(sorted(lst, key=keyfunc), key=keyfunc)]

print desired_list
