def split_line(fn):

    with open(fn, "r") as f:

    	#error because of weird interpretation of files, change later.
    	#input sequence can't have line breaks, 1 complete string
        text = f.readlines() [-1]

    if "\n" in "r":
	   oneliner = fn.replace("\n", " ")
	   print onliner

    # split the text
    words = text.split('-')
    a = 0

    if text.startswith('-'):
        a -= 1

    # for each word in the line:
    for word in words:

        # prints each word on a line
        if word != '':

            print(str(a + 2) + " " + word)
            a = a + len(word) + 1

        else: 
            a = a + 1

#import files w/in same dir pls,
if __name__ == '__main__':
    fn = raw_input('enter file name\n')
    split_line(fn)
