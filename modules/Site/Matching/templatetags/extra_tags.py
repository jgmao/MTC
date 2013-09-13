from django import template
register = template.Library()
@register.filter
def getitem(dict, key): #there are maps inside list, this filter get the idx-th map , and return the value of key    
    return dict[key]

@register.filter
def valueof(list,idx): #there are several maps inside list
    return list[int(idx)].items()

@register.simple_tag
def get_map_val_from_list(list,idx,key):
    return list[int(idx)].get(key)