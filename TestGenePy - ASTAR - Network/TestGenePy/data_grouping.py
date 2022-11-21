def get_groups(data, col):
    return data.groupby(col)

def get_seq_from_group(groups, group_num, index_in_group):
    if group_num in groups.groups.keys():
        group = get_group(groups, group_num)
        return group.iloc[index_in_group, 2]
    else:
        print("Group not found")
        return None

def get_group(groups, group_num):
    if group_num in groups.groups.keys():
        return groups.get_group(group_num)
    else:
        print("Group not found")
        return None

def print_groups(groups):
    for name_of_the_group, group in groups:
        print ("Group = " + str(int(name_of_the_group)))
        print (group)
