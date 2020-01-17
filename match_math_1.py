import sys

query_value=float(sys.argv[1])
upper_query=query_value+(query_value*0.9)
lower_query=query_value-(query_value*0.9)
match_value=float(sys.argv[2])

if match_value < upper_query and match_value > lower_query:
	print(match_value)
