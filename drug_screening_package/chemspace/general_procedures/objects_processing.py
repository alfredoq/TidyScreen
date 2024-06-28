def split_list_in_chunks(list,nbr_items):
    """"
    This function will receive a list of items, and will return a list of lists each one containing the corresponding 'nbr_items'
    ------
    Parameters
    ------
    -list: input list to be processed in chunks
    -nbr_item: number of chunks included in each component list
    """
    list_of_chunks = [list[i:i+nbr_items] for i in range(0, len(list), nbr_items)]

    return list_of_chunks