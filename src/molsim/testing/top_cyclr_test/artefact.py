
<<<<<<< HEAD




#
#
#    while i < len(l):
#        prev_line_empty = False
#        if 'prologue' not in par_order_in:
#            if i < len(l) and not l[i].strip(' ').startswith('['):
#                par_order_in.append('prologue')
#                paragraphs['prologue'] = []
#            while i < len(l) and not l[i].strip(' ').startswith('['):
#                paragraphs['prologue'].append(l[i])
#                i += 1
#            prologue_end = i
#        i += 1
#        # Normal blocks
#        if i < len(l) and l[i].strip(' ').startswith('['):
#            block = l[i].rstrip('\n').strip(' ')
#            if block in par_order_in:
#                block += '_'  # Allows for blocks with the same name
#            par_order_in.append(block)
#            paragraphs[block] = []
#            while i < len(l) and l[i].rstrip('\n').strip(' ') != '':
#                paragraphs[block].append(l[i])
#                i += 1
#            i += 1
#        # # and ; blocks
#        elif i < len(l) and (l[i].strip(' ').startswith('#') or l[i].strip(' ').startswith(';')):
#            i_start = i
#            while i < len(l):
#                if l[i].strip(' ').startswith('#ifdef') or l[i].strip(' ').startswith('#include'):
#                    block = l[i].rstrip('\n').strip(' ')
#                    i = i_start
#                    break
#                elif l[i].rstrip('\n').strip(' ') == '' or i + 1 == len(l):
#                    block = l[i_start].rstrip('\n').strip(' ')
#                    i = i_start
#                    break
#                i += 1
#            par_order_in.append(block)
#            paragraphs[block] = []
#            while i < len(l) and l[i].rstrip('\n').strip(' ') != '':
#                paragraphs[block].append(l[i])
#                i += 1
=======
>>>>>>> 9181c7b (trying to make top_cyclr work with a for loop for line reading.)
