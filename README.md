# Dynamic-Attribute-based-Group-Signature

ABOUT:
- Guide:  Dr. Syed Taqi Ali Sir
- Author: Devashish Arvind Ghate
- Paper Implemented: Dynamic Attribute Based Group Signature with Attribute Anonymity and Tracing in the Standard Model
    by Syed Taqi Ali & B. B. Amberker (https://link.springer.com/chapter/10.1007/978-3-642-41224-0_11)

    Short SUMMARY of the Paper:
      - The aim of the scheme is to allow any member of a certain group to sign a message on behalf of the group,
        but the signer remains anonymous within the group.
      - However, in certain situations, an authority should have the ability to evoke the anonymity of a signer and
        trace the signature.
      - Use in Anonymous attestation, which has practical applications such as in building Trusted Platform Modules
        (TPMs).
      - Attribute Anonymity: Ensures that the attributes used for signing remain hidden, enhancing the privacy of the
        signer.
      - Attribute Tracing: Allows the tracing of attributes used in a signature without revealing the identity of the
        signer.
      - Constant Signature Size: The scheme maintains a fixed signature size, independent of the number of attributes.

    CONCEPTS explained:
    - Pairing based cryptography
        * Cryptography with use of a pairing between elements of two cryptographic groups to a third group with a
        mapping e: G1*G2-> GT.
        * If G1=G2, then the pairing is called SYMMETRIC PAIRING, which is used in this program
        * If G1 != G2, then the pairing is called ASYMMETRIC PAIRING
    - BILINEAR MAP:
        * The map from two groups G1, G2 to a third group GT is the bilinear map.
        (In pbc library, the G1 and G2 groups associated with pairings, are groups of points on an ellitpic curve; GT
        group is currently implemented as a subgroup of a finite field)
        * Denoted (as observed at many places, including the research paper being implemented) with e
    - Properties of bilinear pairing:
        * Bilinearity: e(g^a, h^b) = e(g^b, h^a) = e(g,h)^(ab), where g,h are generators of group G1 & G2 respectively.
        * Non-degenerate: If g,h are generator of G,H then e(g,h) is generator of GT

    The comments/description should be read according to the following:
        - Jump to the main function and start reading the comments
        - Whenever a function call is made, jump to the function body and read comments written just above & inside the
        function body
        - References are mentioned in the format [REF-<Number>]

    IMPORTANT NOTE : Whenever a FUNCTION/DATA-TYPE is encountered for the first time according to the flow mentioned
    above, it is explained in the 'REF.txt' file. For the next time, the Duplicate explanation is avoided. If an
    explanation to FUNCTION/DATA-TYPE is required next time it is being read, please search it in 'REF.txt' file, note
    the reference number & search the explanation of function/data-type by searching [REF-<Number>].

    (It took more effort to avoid duplicate explanation, than it would have taken by writing explanation again.
    Since, it was a good practice to avoid duplicate explanation, I have used such a technique.)

    Following is the WORKING CODE (with comments & descriptions at suitable places) in C programming language for the
    Section "4 Construction" of the paper. The output for the code is present in the file 'output.txt'.

    The code can be compiled using - 'gcc dev.c -lpbc -lgmp' and can be executed using './a.out'
