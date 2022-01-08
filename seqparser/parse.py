import io
from typing import Tuple, Union


class Parser:
    """
    Base Class for Parsing Algorithms
    """
    def __init__(self, filename: str):
        """
        Initialization to be shared by all inherited classes
        
        # What does the `__init__` method do?
            This method will be called immediately upon creating an object. It's a useful
            way to assign baseline attributes of the class (in this case making the filename
            accessible by all methods) but also to run preliminary code or assertions (like 
            checking to see if the file exists at all!)
           
        # Should I ever call the `__init__` method? 
            Like most hidden methods (the double underscored names) this is not generally something
            you call from the outside of a class. However, if you need to specify a different
            `__init__` method for a subclass you will need to call this with the `Super` keyword. We
            won't get into this now, but if you are interested feel free to reach out to the TAs or
            check out the documentation on the `Super` keyword. 
        """
        self.filename = filename

    def get_record(self, f_obj: io.TextIOWrapper) -> Union[Tuple[str, str], Tuple[str, str, str]]:
        """
        Returns a sequencing record that will either be a tuple of two strings (header, sequence)
        or a tuple of three strings (header, sequence, quality). 

        # What's the deal with calling a method by almost the same name?
            it is common in python to see a public method calling a hidden method
            with a similar name. Both of these are accessible to a user (nothing is truly hidden in python)
            but it is a useful way to separate out Class and SubClass specific behavior. 
        
            In this case, we know that the function will return either a tuple of 2 or a tuple of 3. 
            But it is up to the subclass method to define what tuple it will return. 
        
        # Do I need to do this with all my classes?
            Absolutely not. But we want to show you some things you will see often when reading python code
            and give an explanation for why certain practices exist in the language. 

        """
        return self._get_record(f_obj)

    def __iter__(self):
        """
        This is an overriding of the Base Class Iterable function. All classes in python
        have this function, but it is not implemented for all classes in python. 

        # Note on the `__iter__` method
            Generally one doesn't call this method directly as `obj.__iter__()`. Instead it
            lets you use the object itself as an iterable. This is really useful in OOP because it
            allows you to represent and use iterable objects very cleanly. You still can call this
            method directly, but it really takes the fun out of python...

            ## How to use the `__iter__` method
            ```
            parser_obj = Parser(filename)
            for record in parser_obj:
              # do something
            ```

        # Why you should care about generators

            The expected behavior of this function is to create a generator which will lazily load
            the next item in its queue. These are very useful for many bioinformatic tools where you
            don't need everything loaded at once and instead are interested in interacting with the 
            stream (i.e. you need every value once and won't need it again after you use it). This saves
            quite a bit of memory, especially when you are working with billions of sequences and don't 
            need to keep all of them in memory. 
        
        # Distinction between generator functions and other functions
        
            instead of returning a value with the keyword `return`
            a generator must return a value with the keyword `yield`.

            This `yield` keyword will not shortcut the loop it is nested in like a return will
            and instead will pause the loop until the object is taken from it. 
        """

        # The proper way to open a file for reading and writing in python3 is to use the `with` / `as` keywords.
        # and keep the I/O within the nested code block. This will save you from some really nasty bugs that
        # sometimes close a file before everything you expect to be written/read is written/read. 
        # 
        # the interpretation of the following code is that for the lifetime of the filebuffer 
        # returned by the `open` function it will be accessible as the variable `f_obj`
        with open(self.filename, "r") as f_obj:
            
            # this loop will break at some point! 
            # but I will leave it up to you to implement the fix! 

            # You will need to look at the `Try` / `Except` keywords in python
            # and implement an exception for the error you will find in
            # the error message you receive. 
            while True:
                try: 
                    rec = self.get_record(f_obj)
                except:
                    break
                for s in rec: # for sequence in record
                        yield s # yield sequence


    def _get_record(self, f_obj: io.TextIOWrapper) -> Union[Tuple[str, str], Tuple[str, str, str]]:
        """
        a method to be overridden by inherited classes.
        """
        raise NotImplementedError(
                """
                This function is not meant to be called by the Parser Class.
                It is expected to be overridden by `FastaParser` and `FastqParser`""")


class FastaParser(Parser):
    """
    Fasta Specific Parsing
    """
    def _get_record(self, f_obj: io.TextIOWrapper) -> Tuple[str, str]:
        """
        returns the next fasta record
        """

        # First we obtain the header at the top
        for line in f_obj:
            if line[0] == ">": # catching carrot
                seq_name = line[1:].rstrip()
                break
            else:
                break # don't need this, since not getting header

        # at this point we've stored the Header as seq_name
        next(f_obj) # skipping the line since we've obtained the first seq_name/header

        # no longer starts at original seq_name
        for line in f_obj:
            if line[0] == '>': 
                yield seq_name, sequence
                sequence = ''
                seq_name = line[1:].rstrip() # anything after the carrot symbol
                continue
            else:
                sequence = line

        yield seq_name, sequence    

class FastqParser(Parser):
    """
    Fastq Specific Parsing
    """

    def _get_record(self, f_obj: io.TextIOWrapper) -> Tuple[str, str, str]:
        """
        returns the next fastq record
        """

        # get the header
        for line in f_obj:
            if line[0] == '@': # '@' instead of '>'
                seq_name = line[1:].rstrip() # after '@'
            else:
                break

        next(f_obj) #ignore the first header since we have stored already

        for line in f_obj:
            first_line = line
            if line[0] == '@':
                yield, seq_name, sequence, quality
                sequence = ''
                quality = ''
                seq_name = line[1:].rstrip() # everything after '@'
                continue
            elif line[0]== '+': # using this character to determine former and next lines
                sequence = first_line
                quality = next(f_obj) # the line that preceeds
                continue
            else:
                continue

        yield seq_name, sequence, quality



