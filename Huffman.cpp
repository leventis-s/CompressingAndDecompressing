#include "Huffman.h"
#include "map.h"
#include "priorityqueue.h"
#include "error.h"
#include "strlib.h"

using namespace std;

/*
 * Deallocates all nodes in a Huffman tree. We've provided this helper function
 * to you since we also use it in our test driver and figured you might want
 * to make use of it.
 */
void deleteTree(EncodingTreeNode* tree) {
    if (tree != nullptr) {
        deleteTree(tree->zero);
        deleteTree(tree->one);
        delete tree;
    }
}

/*
 * Constructs a Huffman coding tree for the given string, using the algorithm
 * described in lecture.
 *
 * If the input string does not contain at least two different characters,
 * this function reports an error.
 */
EncodingTreeNode* huffmanTreeFor(const string& str) {
   /*Create a map with each char and frequency*/
    Map<char, int> charFrequencies;
    for (char ch: str) {
        if (!charFrequencies.containsKey(ch)) {
            charFrequencies[ch];
        }
        charFrequencies[ch]++;
    }
    if (charFrequencies.size() < 2) {
        error("Less than 2 different characters in string.");
    }
    /*Enqueue all elements as an EncodingTreeNode and their frequency as weight*/
    PriorityQueue<EncodingTreeNode*> characters;
    for (char key: charFrequencies) {
        EncodingTreeNode* element = new EncodingTreeNode;
        element->ch = key;
        element->zero = nullptr;
        element->one = nullptr;
        characters.enqueue(element, charFrequencies[key]);
    }
    /* Dequeue 2 elemnts of least frequency and then enqueue them as one larger node*/
    while (characters.size() > 1) {
        EncodingTreeNode* parent = new EncodingTreeNode;
        int weightZero = characters.peekPriority();
        parent->zero = characters.dequeue();
        int weightOne = characters.peekPriority();
        parent->one = characters.dequeue();
        characters.enqueue(parent, weightZero + weightOne);
    }
    return characters.dequeue(); //returns pointer to final node
}

/*
 * Given a Queue<Bit> containing a compressed message and a tree that was used
 * to encode those bits, decodes the bits back to the original message.
 */
string decodeText(Queue<Bit>& bits, EncodingTreeNode* tree) {
    Bit b = 0;
    string result = "";
    while (!bits.isEmpty()) {
        EncodingTreeNode* decoder = tree;
        while (decoder->zero != nullptr && decoder->one != nullptr) { //While it's not leaf node
            b = bits.dequeue();
            if (b == 0) {
                decoder = decoder->zero;
            }
            else {
                decoder = decoder->one;
            }
        }
        result += decoder->ch;
    }
    return result;
}

/*
 * Given the tree, it creates a map of each char and its code.
 */
void createTable(EncodingTreeNode* tree, Map<char, string>& table, string soFar) {
    if (tree->zero == nullptr) {
        table[tree->ch] = soFar;
    }
    else {
        createTable(tree->zero, table, soFar + "0");
        createTable(tree->one, table, soFar + "1");
    }
}

/*
 * Given a string and a Huffman encoding tree, encodes that text using the tree
 * and outputs a Queue<Bit> corresponding to the encoded representation.
 */
Queue<Bit> encodeText(const string& str, EncodingTreeNode* tree) {
    Queue<Bit> result;
    Map<char, string> table;
    createTable(tree, table, "");
    /* Decode string using a map */
    for (char ch: str) {
        string code = table[ch];
        for (char num: code) {
            result.enqueue(charToInteger(num));
        }
    }
    return result;
}

/*
 * Recursively creates the tree using a queue of its bits and a seperate queue of its leaves.
 */
void decodeTreeRec(Queue<Bit>& bits, Queue<char>& leaves, EncodingTreeNode* decodeTree) {
    if (bits.isEmpty()) {
        return;
    }
    Bit b = bits.dequeue();
    /* B determines if it is a leaf or not */
    if (b == 0) {
        decodeTree->ch = leaves.dequeue();
        decodeTree->zero = nullptr;
        decodeTree->one = nullptr;
    }
    else {
        EncodingTreeNode* zeroNode = new EncodingTreeNode;
        decodeTree->zero = zeroNode;
        decodeTreeRec(bits, leaves, zeroNode);
        EncodingTreeNode* oneNode = new EncodingTreeNode;
        decodeTree->one = oneNode;
        decodeTreeRec(bits, leaves, oneNode);
    }
}

/*
 * Decodes the given Queue<Bit> and Queue<char> into a Huffman coding tree.
 */
EncodingTreeNode* decodeTree(Queue<Bit>& bits, Queue<char>& leaves) {
    EncodingTreeNode* decodeTree = new EncodingTreeNode;
    decodeTreeRec(bits, leaves, decodeTree);
    return decodeTree;
}


/*
 * Encodes the given Huffman tree as a Queue<Bit> and Queue<char> in the manner
 * specified in the assignment handout.
 */
void encodeTree(EncodingTreeNode* tree, Queue<Bit>& bits, Queue<char>& leaves) {
    /* If it has a nullptr it is a leaf node */
    if (tree->zero == nullptr) {
        bits.enqueue(0);
        leaves.enqueue(tree->ch);
    }
    else {
        bits.enqueue(1);
        encodeTree(tree->zero, bits, leaves);
        encodeTree(tree->one, bits, leaves);
    }
}

/*
 * Compresses the given text string using Huffman coding, producing as output
 * a HuffmanResult containing the encoded tree and message.
 */
HuffmanResult compress(const string& text) {
    HuffmanResult result;
    EncodingTreeNode* tree = huffmanTreeFor(text);
    result.messageBits = encodeText(text, tree);
    Queue<Bit> bits;
    Queue<char> leaves;
    encodeTree(tree, bits, leaves);
    result.treeLeaves = leaves;
    result.treeBits = bits;
    deleteTree(tree);
    return result;
}

/*
 * Decompresses the given HuffmanResult and returns the string it represents.
 */
string decompress(HuffmanResult& file) {
    EncodingTreeNode* tree = decodeTree(file.treeBits, file.treeLeaves);
    string result = decodeText(file.messageBits, tree);
    deleteTree(tree);
    return result;
}

/* * * * * * Test Cases Below This Point * * * * * */
#include "GUI/SimpleTest.h"

/* TODO: Add your own custom tests here! */

bool isEncodingTree(EncodingTreeNode* tree) {
    /* The empty tree is not a Huffman tree. */
    if (tree == nullptr) return false;

    /* If we have one missing child, we should have two missing children. */
    if ((tree->zero == nullptr) != (tree->one == nullptr)) return false;

    /* If we have children at all, they need to be Huffman trees. */
    return tree->zero == nullptr || (isEncodingTree(tree->zero) && isEncodingTree(tree->one));
}

/* Utility function to test if two trees are equal. This is adapted from Section
 * Handout 8 and particularized to Huffman trees.
 */
bool areEqual(EncodingTreeNode* lhs, EncodingTreeNode* rhs) {
    /* Base case: If either is a leaf, both should be. */
    bool lhsLeaf = lhs->zero == nullptr && lhs->one == nullptr;
    bool rhsLeaf = rhs->zero == nullptr && rhs->one == nullptr;
    if (lhsLeaf || rhsLeaf) {
        return lhs->ch == rhs->ch && lhsLeaf == rhsLeaf;
    }

    /* Otherwise, they're both internal nodes. Check that they match. */
    return areEqual(lhs->zero, rhs->zero) && areEqual(lhs->one, rhs->one);
}

/* Utility function to return a string of all possible characters. */
string pangrammaticString() {
    string result;

    char ch = numeric_limits<char>::min();
    result += ch;
    do {
        ch++;
        result += ch;
    } while (ch != numeric_limits<char>::max());

    return result;
}

/* Utility function that makes an inefficient but still valid encoding tree
 * for the given characters.
 */

EncodingTreeNode* strandTreeFor(const string& text, size_t index = 0) {
    if (index == text.size()) error("No characters provided to strandTreeFor.");

    /* We always get a leaf node. */
    EncodingTreeNode* leaf = new EncodingTreeNode {
        text[index], nullptr, nullptr
    };

    /* Last character? If so, that's all. */
    if (index + 1 == text.size()) return leaf;

    /* Otherwise, build a larger tree. */
    else return new EncodingTreeNode {
        ' ', leaf, strandTreeFor(text, index + 1)
    };
}

STUDENT_TEST("huffmanTreeFor builds tree for same characters that are non-consecutive.") {
    EncodingTreeNode* reference = new EncodingTreeNode {
        ' ', new EncodingTreeNode {'a', nullptr, nullptr}, new EncodingTreeNode {'b', nullptr, nullptr}
    };

    EncodingTreeNode* tree = huffmanTreeFor("abbaabb");
    EXPECT(isEncodingTree(tree));
    EXPECT(areEqual(tree, reference));

    deleteTree(reference);
    deleteTree(tree);
}

STUDENT_TEST("huffmanTreeFor builds tree for characters with same weight.") {
    EncodingTreeNode* reference = new EncodingTreeNode {
        ' ', new EncodingTreeNode {'b', nullptr, nullptr}, new EncodingTreeNode {'a', nullptr, nullptr}
    };

    EncodingTreeNode* tree = huffmanTreeFor("abbaab");
    EXPECT(isEncodingTree(tree));
    EXPECT(areEqual(tree, reference));

    deleteTree(reference);
    deleteTree(tree);
}

STUDENT_TEST("encodeText works on sample of just two leaf nodes.") {
    /* This tree:
     *                 *
     *                / \
     *               O   N
     */
    EncodingTreeNode* tree = new EncodingTreeNode {
           '*',
           new EncodingTreeNode { 'O', nullptr, nullptr },
           new EncodingTreeNode{ 'N', nullptr, nullptr }
    };

    /* What you get if you encode NOON with this tree. */
    Queue<Bit> expected = { 1, 0, 0, 1};

    EXPECT_EQUAL(encodeText("NOON", tree), expected);

    deleteTree(tree);
}

STUDENT_TEST("decodeText works on tree with just 2 nodes.") {
    /* This tree:
     *                 *
     *                / \
     *               O   N
     */
    EncodingTreeNode* tree = new EncodingTreeNode {
        '*',
            new EncodingTreeNode { 'O', nullptr, nullptr },
            new EncodingTreeNode { 'N',
                new EncodingTreeNode{ '*',
            },
        },
    };

    Queue<Bit> bits = { 0, 1, 1, 0};

    EXPECT_EQUAL(decodeText(bits, tree), "ONNO");

    deleteTree(tree);
}

STUDENT_TEST("Can decode a large tree.") {
    /* This encodes this tree:
     *
     *                 *
     *               /   \
     *              *     *
     *             / \   / \
     *            A   B C   *
     *                     / \
     *                    D   E
     */
    Queue<Bit>  bits   = { 1, 1, 0, 0, 1, 0, 1, 0, 0};
    Queue<char> leaves = { 'A', 'B', 'C', 'D', 'E'};

    EncodingTreeNode* tree = decodeTree(bits, leaves);
    EXPECT(isEncodingTree(tree));

    /* Confirm this is the right tree. */
    EncodingTreeNode* expected = new EncodingTreeNode {
        '*',
            new EncodingTreeNode {
                '*',
                new EncodingTreeNode { 'A', nullptr, nullptr },
                new EncodingTreeNode { 'B', nullptr, nullptr },
            },
            new EncodingTreeNode {
                '*',
                new EncodingTreeNode { 'C', nullptr, nullptr },
                new EncodingTreeNode {
                    '*',
                    new EncodingTreeNode { 'D', nullptr, nullptr },
                    new EncodingTreeNode { 'E', nullptr, nullptr },
                },
            },
    };

    EXPECT(areEqual(tree, expected));

    deleteTree(tree);
    deleteTree(expected);
}

STUDENT_TEST("Can encode a tree with just 2 leaf nodes of equal weight.") {
    /* Build an encoding tree for "AC" It should look like this:
     *
     *                 *
     *                / \
     *               A   C
     *
     * This will compress down to
     *
     *        100
     *        AC
     */
    EncodingTreeNode* tree = huffmanTreeFor("AC");

    Queue<Bit>  bits;
    Queue<char> leaves;

    encodeTree(tree, bits, leaves);

    Queue<Bit>  expectedBits   = { 1, 0, 0};
    Queue<char> expectedLeaves = { 'C', 'A' };

    EXPECT_EQUAL(bits,   expectedBits);
    EXPECT_EQUAL(leaves, expectedLeaves);

    deleteTree(tree);
}

STUDENT_TEST("Can decompress a file with minimum size tree.") {
    HuffmanResult file = {
        { 1, 0, 0 },
        { 'a', 'c'},
        { 0, 1 }
    };

    EXPECT_EQUAL(decompress(file), "ac");
}

/* * * * * Provided Tests Below This Point * * * * */
#include <limits>

/* Utility function to test if a purported Huffman tree is indeed a Huffman tree.
 * Specifically, this checks that each internal node has either zero or two
 * children. There are other ways you could produce an invalid Huffman tree - for
 * example, by having uninitialized pointers or by linking in a cycle - but we
 * don't test for that here.
 */


PROVIDED_TEST("huffmanTreeFor reports errors on invalid inputs.") {
    EXPECT_ERROR(huffmanTreeFor(""));    // No characters
    EXPECT_ERROR(huffmanTreeFor("a"));   // One character
    EXPECT_ERROR(huffmanTreeFor("aaa")); // One character
}

PROVIDED_TEST("huffmanTreeFor builds tree for two characters.") {
    EncodingTreeNode* reference = new EncodingTreeNode {
        ' ', new EncodingTreeNode {'a', nullptr, nullptr}, new EncodingTreeNode {'b', nullptr, nullptr}
    };

    EncodingTreeNode* tree = huffmanTreeFor("aaabbbb");
    EXPECT(isEncodingTree(tree));
    EXPECT(areEqual(tree, reference));

    deleteTree(reference);
    deleteTree(tree);
}

PROVIDED_TEST("huffmanTreeFor works on the full range of characters.") {
    /* Get a string of all possible characters, then pair them off and see what we find. */
    string allChars = pangrammaticString();
    for (size_t i = 0; i < allChars.size(); i += 2) {
        string toEncode;
        toEncode += allChars[i];
        toEncode += allChars[i + 1];
        toEncode += allChars[i + 1];

        EncodingTreeNode* reference = new EncodingTreeNode {
            ' ',
            new EncodingTreeNode {allChars[i], nullptr, nullptr},
            new EncodingTreeNode {allChars[i + 1], nullptr, nullptr}
        };

        EncodingTreeNode* tree = huffmanTreeFor(toEncode);
        EXPECT(isEncodingTree(tree));
        EXPECT(areEqual(tree, reference));

        deleteTree(reference);
        deleteTree(tree);
    }
}

PROVIDED_TEST("huffmanTreeFor uses cumulative weights (v1).") {
    /* This tree:
     *                 *
     *                / \
     *               *   D
     *              / \
     *             C   *
     *                / \
     *               A   B
     */
    EncodingTreeNode* reference = new EncodingTreeNode {
        '*',
            new EncodingTreeNode { '!',
                new EncodingTreeNode { 'C', nullptr, nullptr },
                new EncodingTreeNode { '?',
                    new EncodingTreeNode { 'A', nullptr, nullptr },
                    new EncodingTreeNode { 'B', nullptr, nullptr }
                }
            },
            new EncodingTreeNode { 'D', nullptr, nullptr }
    };

    /* Ax2, Bx3, Cx4, Dx10 */
    EncodingTreeNode* tree = huffmanTreeFor("AABBBCCCCDDDDDDDDDD");
    EXPECT(isEncodingTree(tree));
    EXPECT(areEqual(tree, reference));

    deleteTree(reference);
    deleteTree(tree);
}

PROVIDED_TEST("huffmanTreeFor uses cumulative weights (v2).") {
    /*
     *          *
     *       /     \
     *      *       *
     *     / \     / \
     *    D   E   F   *
     *               / \
     *              C   *
     *                 / \
     *                A   B
     */
    EncodingTreeNode* reference =new EncodingTreeNode {
        ' ',
        new EncodingTreeNode {
            ' ',
            new EncodingTreeNode{ 'D', nullptr, nullptr },
            new EncodingTreeNode{ 'E', nullptr, nullptr }
        },
        new EncodingTreeNode {
            ' ',
            new EncodingTreeNode{ 'F', nullptr, nullptr },
            new EncodingTreeNode {
                ' ',
                new EncodingTreeNode{ 'C', nullptr, nullptr },
                new EncodingTreeNode{
                    ' ',
                    new EncodingTreeNode{ 'A', nullptr, nullptr },
                    new EncodingTreeNode{ 'B', nullptr, nullptr },
                }
            }
        }
    };

    /* Ax2, Bx3, Cx4, Dx6, Ex7, Fx8 */
    EncodingTreeNode* tree = huffmanTreeFor("AABBBCCCCDDDDDDEEEEEEEFFFFFFFF");
    EXPECT(isEncodingTree(tree));

    EXPECT(areEqual(tree, reference));

    deleteTree(reference);
    deleteTree(tree);
}

PROVIDED_TEST("decodeText works on small sample.") {
    /* This tree:
     *                 *
     *                / \
     *               O   *
     *                  / \
     *                 *   N
     *                / \
     *               M   S
     */
    EncodingTreeNode* tree = new EncodingTreeNode {
        '*',
            new EncodingTreeNode { 'O', nullptr, nullptr },
            new EncodingTreeNode { '*',
                new EncodingTreeNode{ '*',
                    new EncodingTreeNode { 'M', nullptr, nullptr },
                    new EncodingTreeNode { 'S', nullptr, nullptr }
                },
                new EncodingTreeNode{ 'N', nullptr, nullptr }
            }
    };

    /* What you get if you encode MONSOON with this tree. */
    Queue<Bit> bits = { 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1 };

    EXPECT_EQUAL(decodeText(bits, tree), "MONSOON");

    deleteTree(tree);
}

PROVIDED_TEST("Can decode all char values.") {
    /* All possible characters. */
    string allChars = pangrammaticString();

    /* Try decoding all pairs of adjacent characters. */
    for (size_t i = 0; i < allChars.size(); i += 2) {
        string expected;
        expected += allChars[i];
        expected += allChars[i + 1];
        expected += allChars[i + 1];

        EncodingTreeNode* tree = new EncodingTreeNode {
            ' ',
            new EncodingTreeNode {allChars[i], nullptr, nullptr},
            new EncodingTreeNode {allChars[i + 1], nullptr, nullptr}
        };

        /* Decode the bitstream 011, which should map back to the expected
         * string.
         */
        Queue<Bit> bits = { 0, 1, 1 };
        EXPECT_EQUAL(decodeText(bits, tree), expected);

        deleteTree(tree);
    }
}

PROVIDED_TEST("encodeText works on small sample.") {
    /* This tree:
     *                 *
     *                / \
     *               O   *
     *                  / \
     *                 *   N
     *                / \
     *               M   S
     */
    EncodingTreeNode* tree = new EncodingTreeNode {
           '*',
           new EncodingTreeNode { 'O', nullptr, nullptr },
               new EncodingTreeNode { '*',
               new EncodingTreeNode{ '*',
               new EncodingTreeNode { 'M', nullptr, nullptr },
               new EncodingTreeNode { 'S', nullptr, nullptr }
            },
            new EncodingTreeNode{ 'N', nullptr, nullptr }
        }
    };

    /* What you get if you encode MONSOON with this tree. */
    Queue<Bit> expected = { 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1 };

    EXPECT_EQUAL(encodeText("MONSOON", tree), expected);

    deleteTree(tree);
}

PROVIDED_TEST("Can encode all char values.") {
    /* All possible characters. */
    string allChars = pangrammaticString();

    /* Try encoding all pairs of adjacent characters. */
    for (size_t i = 0; i < allChars.size(); i += 2) {
        string toEncode;
        toEncode += allChars[i];
        toEncode += allChars[i + 1];
        toEncode += allChars[i + 1];

        EncodingTreeNode* tree = new EncodingTreeNode {
                ' ',
                new EncodingTreeNode {allChars[i], nullptr, nullptr},
                new EncodingTreeNode {allChars[i + 1], nullptr, nullptr}
        };

        /* See what bits we get back. We should get 011, since the first
         * character has code 0 and the second has code 1.
         */
        Queue<Bit> bits = encodeText(toEncode, tree);
        Queue<Bit> expected = { 0, 1, 1 };

        EXPECT_EQUAL(bits, expected);

        deleteTree(tree);
    }
}

PROVIDED_TEST("decodeText undoes encodeText on range of sample strings.") {
    Vector<string> testCases = {
        "THAT THAT IS IS THAT THAT IS NOT IS NOT IS THAT IT IT IS",
        "AABAAABBABAAABAAAA",
        ":-) :-D XD <(^_^)>",
        pangrammaticString(),
    };

    for (string test: testCases) {
        /* Use a silly encoding scheme, but one that works regardless of what
         * characters are provided.
         */
        EncodingTreeNode* tree = strandTreeFor(test);
        EXPECT(isEncodingTree(tree));

        Queue<Bit> bits = encodeText(test, tree);
        string result = decodeText(bits, tree);

        EXPECT_EQUAL(test.size(), result.size());
        EXPECT_EQUAL(test, result);

        deleteTree(tree);
    }
}

PROVIDED_TEST("Can decode an example tree.") {
    /* This encodes this tree:
     *
     *                 *
     *                / \
     *               *   C
     *              / \
     *             A   B
     */
    Queue<Bit>  bits   = { 1, 1, 0, 0, 0 };
    Queue<char> leaves = { 'A', 'B', 'C' };

    EncodingTreeNode* tree = decodeTree(bits, leaves);
    EXPECT(isEncodingTree(tree));

    /* Confirm this is the right tree. */
    EncodingTreeNode* expected = new EncodingTreeNode {
        '*',
            new EncodingTreeNode {
                '*',
                new EncodingTreeNode { 'A', nullptr, nullptr },
                new EncodingTreeNode { 'B', nullptr, nullptr },
            },
            new EncodingTreeNode { 'C', nullptr, nullptr }
    };

    EXPECT(areEqual(tree, expected));

    deleteTree(tree);
    deleteTree(expected);
}

PROVIDED_TEST("Can decode trees using all possible char values.") {
    /* All possible characters. */
    string allChars = pangrammaticString();

    /* Try encoding all pairs of adjacent characters. */
    for (size_t i = 0; i < allChars.size(); i += 2) {
        EncodingTreeNode* expected = new EncodingTreeNode {
            ' ',
            new EncodingTreeNode {allChars[i], nullptr, nullptr},
            new EncodingTreeNode {allChars[i + 1], nullptr, nullptr}
        };
        Queue<Bit>  treeBits   = { 1, 0, 0 };
        Queue<char> treeLeaves = { allChars[i], allChars[i + 1] };

        EncodingTreeNode* tree = decodeTree(treeBits, treeLeaves);
        EXPECT(isEncodingTree(tree));
        EXPECT(areEqual(tree, expected));

        deleteTree(tree);
        deleteTree(expected);
    }
}

PROVIDED_TEST("Can encode an example tree.") {
    /* Build an encoding tree for "ABBCCCC." It should look like this:
     *
     *                 *
     *                / \
     *               *   C
     *              / \
     *             A   B
     *
     * This will compress down to
     *
     *        11000
     *        ABC
     */
    EncodingTreeNode* tree = huffmanTreeFor("ABBCCCC");

    Queue<Bit>  bits;
    Queue<char> leaves;

    encodeTree(tree, bits, leaves);

    Queue<Bit>  expectedBits   = { 1, 1, 0, 0, 0 };
    Queue<char> expectedLeaves = { 'A', 'B', 'C' };

    EXPECT_EQUAL(bits,   expectedBits);
    EXPECT_EQUAL(leaves, expectedLeaves);

    deleteTree(tree);
}

PROVIDED_TEST("Can encode trees using all possible char values.") {
    /* All possible characters. */
    string allChars = pangrammaticString();

    /* Try encoding all pairs of adjacent characters. */
    for (size_t i = 0; i < allChars.size(); i += 2) {
        EncodingTreeNode* tree = new EncodingTreeNode {
            ' ',
            new EncodingTreeNode {allChars[i], nullptr, nullptr},
            new EncodingTreeNode {allChars[i + 1], nullptr, nullptr}
        };

        /* See what we get back. It should be the bitstring 100 (root with
         * two children) and the two leaves, in order.
         */
        Queue<Bit>  treeBits;
        Queue<char> treeLeaves;

        Queue<Bit>  expectedBits = { 1, 0, 0 };
        Queue<char> expectedLeaves = { allChars[i], allChars[i + 1] };

        encodeTree(tree, treeBits, treeLeaves);
        EXPECT_EQUAL(treeBits, expectedBits);
        EXPECT_EQUAL(treeLeaves, expectedLeaves);

        deleteTree(tree);
    }
}

PROVIDED_TEST("decodeTree undoes encodeTree on sample strings.") {
    /* Make a Huffman tree for the string of all characters. */
    EncodingTreeNode* sourceTree = huffmanTreeFor(pangrammaticString());
    EXPECT(isEncodingTree(sourceTree));

    /* Encode, then decode it. */
    Queue<Bit>  bits;
    Queue<char> leaves;
    encodeTree(sourceTree, bits, leaves);

    EncodingTreeNode* resultTree = decodeTree(bits, leaves);
    EXPECT(isEncodingTree(resultTree));
    EXPECT(areEqual(sourceTree, resultTree));

    deleteTree(sourceTree);
    deleteTree(resultTree);
}

PROVIDED_TEST("Can decompress a small sample file.") {
    HuffmanResult file = {
        { 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0 },
        { 'u', 'k', 'p', 'n', 'a', 'm', 'h' },
        { 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1,
          0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0,
          0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0 }
    };

    EXPECT_EQUAL(decompress(file), "humuhumunukunukuapuaa");
}

PROVIDED_TEST("Compress reports errors on bad inputs.") {
    EXPECT_ERROR(compress(""));
    EXPECT_ERROR(compress("A"));
    EXPECT_ERROR(compress("AAAA"));
}

PROVIDED_TEST("Can compress a small sample file.") {
    HuffmanResult file = compress("ABANANAABANDANA");
    Queue<Bit>  treeBits    = { 1, 1, 1, 0, 0, 0, 0 };
    Queue<char> treeChars   = { 'D', 'B', 'N', 'A' };
    Queue<Bit>  messageBits = { 1, 0, 0, 1, 1, 0, 1, 1, 0,
                                1, 1, 1, 0, 0, 1, 1, 0, 1,
                                0, 0, 0, 1, 0, 1, 1 };

    EXPECT_EQUAL(file.treeBits, treeBits);
    EXPECT_EQUAL(file.treeLeaves, treeChars);
    EXPECT_EQUAL(file.messageBits, messageBits);
}

PROVIDED_TEST("Compress undoes decompress on a range of strings.") {
    Vector<string> testCases = {
        "THAT THAT IS IS THAT THAT IS NOT IS NOT IS THAT IT IT IS",
        "AABAAABBABAAABAAAA",
        ":-) :-D XD <(^_^)>",
        pangrammaticString(),
    };

    for (string test: testCases) {
        HuffmanResult file = compress(test);
        string result = decompress(file);

        EXPECT_EQUAL(result.size(), test.size());
        EXPECT_EQUAL(test, result);
    }
}
