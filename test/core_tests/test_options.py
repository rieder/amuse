from amuse.support.options import OptionalAttributes, options, option, GlobalOptions
from amuse.test import amusetest
import StringIO
import textwrap
import os

class OptionsTestsClass(OptionalAttributes):
    option_sections = ('amuse',)
    
    def string_option(self):
        return "default string"
    
    def boolean_option(self):
        return False
        
    def int_option(self):
        return 2
        
    def float_option(self):
        return 2.5
        
    
class OptionsTests(amusetest.TestCase):
    ini_contents = textwrap.dedent("""
    [amuse]
    string_option=a string
    boolean_option=1
    int_option=1
    float_option=1.5
    [bhtree]
    mode=gpu
    int_option=3
    """)
    
    
    def test1(self):
        global_options = GlobalOptions()
        global_options.config.readfp(StringIO.StringIO(self.ini_contents))
        instance = OptionsTestsClass()
        
        testoption = option(OptionsTestsClass.string_option, global_options=global_options)
        self.assertEquals(testoption.__get__(instance, type(instance)), "a string")
        
    def test2(self):
        global_options = GlobalOptions()
        global_options.config.readfp(StringIO.StringIO(self.ini_contents))
        instance = OptionsTestsClass()
        
        testoption = option(OptionsTestsClass.int_option, type="int", global_options=global_options)
        self.assertEquals(testoption.__get__(instance, type(instance)), 1)
        
        testoption = option(OptionsTestsClass.boolean_option, type="boolean", global_options=global_options)
        self.assertEquals(testoption.__get__(instance, type(instance)), True)
        
        testoption = option(OptionsTestsClass.float_option, type="float", global_options=global_options)
        self.assertEquals(testoption.__get__(instance, type(instance)), 1.5) 
        
    def test3(self):
        global_options = GlobalOptions()
        global_options.config.readfp(StringIO.StringIO(self.ini_contents))
        instance = OptionsTestsClass()
        instance.option_sections =('unknown')
        
        testoption = option(OptionsTestsClass.string_option, global_options=global_options)
        self.assertEquals(testoption.__get__(instance, type(instance)), "default string")
        
        testoption = option(OptionsTestsClass.int_option, type="int", global_options=global_options)
        self.assertEquals(testoption.__get__(instance, type(instance)), 2)
        
        testoption = option(OptionsTestsClass.boolean_option, type="boolean", global_options=global_options)
        self.assertEquals(testoption.__get__(instance, type(instance)), False)
        
        testoption = option(OptionsTestsClass.float_option, type="float", global_options=global_options)
        self.assertEquals(testoption.__get__(instance, type(instance)), 2.5)
        
    
    def test4(self):
        global_options = GlobalOptions()
        global_options.config.readfp(StringIO.StringIO(self.ini_contents))
        instance = OptionsTestsClass()
        
        testoption = option(OptionsTestsClass.int_option, type="int", sections=("bhtree",), global_options=global_options)
        self.assertEquals(testoption.__get__(instance, type(instance)), 1)
        
        instance.option_sections =('unknown')
        self.assertEquals(testoption.__get__(instance, type(instance)), 3)
        
        
    
    def test5(self):
        global_options = GlobalOptions()
        global_options.config.readfp(StringIO.StringIO(self.ini_contents))
        instance = OptionsTestsClass()
        instance.option_sections =('unknown')
        
        testoption = option(type="int", global_options=global_options)
        testoption = testoption(OptionsTestsClass.int_option)
        
        self.assertEquals(testoption.__get__(instance, type(instance)), 2)
        
    def test6(self):
        global_options = GlobalOptions()
        global_options.config.readfp(StringIO.StringIO(self.ini_contents))
        class DecoratedMethods(OptionalAttributes):
            option_sections =('unknown')
            
            @option(type="int", global_options=global_options)
            def int_option(self):
                return self.i
            
                
        instance = DecoratedMethods()
        instance.i = 10
        self.assertEquals(instance.int_option, 10)
        
    def test7(self):
        global_options = GlobalOptions()
        global_options.config.readfp(StringIO.StringIO(self.ini_contents))
        class DecoratedMethods(OptionalAttributes):
            option_sections =('unknown')
            i = 8
            def __init__(self, **optional_arguments):
                OptionalAttributes.__init__(self, **optional_arguments)
                
            @option(type="int", global_options=global_options)
            def int_option(self):
                return self.i
            
                
        instance = DecoratedMethods()
        instance.i = 10
        instance.int_option = 12
        self.assertEquals(instance.int_option, 12)
        
        instance = DecoratedMethods(int_option=14)
        self.assertEquals(instance.int_option, 14)
        

    def test8(self):
        global_options = GlobalOptions()
        global_options.config.readfp(StringIO.StringIO(self.ini_contents))
        class A(OptionalAttributes):
            option_sections = ('amuse')
            
            def __init__(self, **optional_arguments):
                OptionalAttributes.__init__(self, **optional_arguments)
            
            @option(type="int", global_options=global_options)
            def int_option_a(self):
                return self.i
                
        class B(OptionalAttributes):
            option_sections =('amuse')
            
            def __init__(self, **optional_arguments):
                OptionalAttributes.__init__(self, **optional_arguments)
        
            @option(type="int", global_options=global_options)
            def int_option_b(self):
                return 10
                
        class C(A,B):
            
            def __init__(self, **optional_arguments):
                A.__init__(self, **optional_arguments)
                B.__init__(self, **optional_arguments)
                
        
        instance = C(int_option_a=14)
        self.assertEquals(instance.int_option_a, 14)
        self.assertEquals(instance.int_option_b, 10)
        
    def test9(self):
        global_options = GlobalOptions()
        global_options.config.readfp(StringIO.StringIO(self.ini_contents))
        class A(OptionalAttributes):
            option_sections = ('amuse')
            
            def __init__(self, **optional_arguments):
                OptionalAttributes.__init__(self, **optional_arguments)
            
            @option(type="int", global_options=global_options)
            def int_option_a(self):
                return self.i
                
        class B(OptionalAttributes):
            option_sections =('amuse')
            
            def __init__(self, **optional_arguments):
                OptionalAttributes.__init__(self, **optional_arguments)
        
            @option(type="int", global_options=global_options)
            def int_option_b(self):
                return 10
                
        class C(A,B):
            
            def __init__(self, **optional_arguments):
                A.__init__(self, **optional_arguments)
                B.__init__(self, **optional_arguments)
                
        
        instance = C(int_option_a=14)
        self.assertEquals(instance.int_option_a, 14)
        self.assertEquals(instance.int_option_b, 10)
        
    def test10(self):
        global_options = GlobalOptions()
        self.assertEquals(global_options.rcfilename, 'amuserc')
        global_options.rcfilename = 'test_for_amuserc'
        with open(global_options.rcfilename, "w") as f:
            f.write(self.ini_contents)
            
        self.assertTrue(global_options.rcfilepath.endswith('test_for_amuserc'))
        global_options.load()
        
        instance = OptionsTestsClass()
        testoption = option(OptionsTestsClass.int_option, type="int", global_options=global_options)
        self.assertEquals(testoption.__get__(instance, type(instance)), 1)
        
        os.remove(global_options.rcfilepath)
        
    
    def test11(self):
        global_options = GlobalOptions()
        global_options.config.readfp(StringIO.StringIO(self.ini_contents))
        class DecoratedMethods(OptionalAttributes):
            option_sections =('amuse',)
           
            def __init__(self, **optional_arguments):
                OptionalAttributes.__init__(self, **optional_arguments)
                
            @option(choices=("1","2"), global_options=global_options)
            def int_option(self):
                return 2
            
                
        instance = DecoratedMethods()
        self.assertEquals(instance.int_option, "1")
        all_options = list(instance.iter_options())
        self.assertEquals(len(all_options), 1)
        self.assertEquals(all_options[0].name, "int_option")

    def test12(self):
        global_options = GlobalOptions()
        global_options.read_from_ini_string(self.ini_contents)
        print global_options.to_ini_string()
        ini_string = global_options.to_ini_string()
        
        self.assertTrue(ini_string.find("= a string") > 0)
    
    
