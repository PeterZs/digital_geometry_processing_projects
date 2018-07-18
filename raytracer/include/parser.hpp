#ifndef PARSER_H
#define PARSER_H

#include "scene.hpp"

#include <climits>
#include <deque>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <memory>

// Scene file parser
// The different types of tokens that can be lexed.
enum TokenType {
	STRING,
	NUMBER,
	NAME,
	ARRAY_BEGIN,
	ARRAY_END,
	END_OF_FILE,
	ERROR
};


// A class to represent a single token that has been lexed.
class Token 
{
public:
	TokenType type;

	// Variables for the different data types.
	double number;
	std::string string;

	// A few constructors for directly assigning values.
	Token(TokenType type) : type(type) {}
	Token(double value) : type(NUMBER), number(value) {}
	Token(TokenType type, std::string value) : type(type), string(value) {}

	// Equality operator.
	bool operator==(Token const &other) const;
};


// Convenience to write a Token directly to an output stream.
std::ostream& operator<<(std::ostream &out, Token const &token);


// The lexing class itself.
class Lexer 
{
private:
	// The input stream to read from.
	std::istream& _input;

	// Raw function to read the input stream and return the next token.
	Token _processStream();

	// A temporary storage for lexed tokens that have not been parsed yet.
	std::deque<Token> _buffer;

public:
	// Construct a lexer. Must be given an input stream.
	Lexer(std::istream& input) : _input(input) {}

	// Peek at the next token, but don't consume it.
	Token peek(unsigned int index = 0);

	// Get the next token.
	Token next();

	// Skip a number of tokens.
	void skip(unsigned int count = 1);

	// The remaining functions will throw a std::string exception if the
	// are unable to perform as requested.

	// Retrieve a command name.
	std::string getName();

	// Get a list of numbers. Min/max refer to the required size of the
	// list.
	std::vector<double> getNumbers(unsigned int min = 0,
		unsigned int max = UINT_MAX);

	// Get a single number.
	double getNumber();

	// Get a single string.
	std::string getString();

	// Get a parameter list (i.e. strings mapping to number lists).
	// Min/max apply to each number list as for getNumbers().
	ParamList getParamList(unsigned int min = 0,
		unsigned int max = UINT_MAX);
};


class Parser 
{
private:
    Lexer lexer; // The lexer we will use to break up the file into tokens.
    
    std::vector<Matrix> transformStack;  // Transformation stack.

    // The following functions parse all of the commands that may be found
    // in a .ray file.
    void parseDimensions();
    void parsePerspective();
    void parseLookAt();
    void parseMaterial();

    void parsePushMatrix();
    void parsePopMatrix();
    void parseTranslate();
    void parseRotate();
    void parseScale();

    void parseSphere();
    void parsePlane();
    void parseMesh();
    void parseConic();

    void parsePointLight();

    // Parses the common parts of each object (trailing material name),
    // and sets up the object in the scene.
    void finishObject(Object *obj);

public:
    Scene scene; // The scene that will be created when parsing.

    Parser(std::istream &input) : lexer(input) {}

    // Parse the file or stream passed into the constructor.
    // Store the results in scene.
    // Returns false on failure; will print to std::cerr.
    bool parse();
};


#endif
