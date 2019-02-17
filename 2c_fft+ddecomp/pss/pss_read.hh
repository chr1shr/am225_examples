#ifndef PSS_READ_HH
#define PSS_READ_HH

#include <vector>

class pss_read { 
    public:
        /** The side length of the perfect squared square. */
        const int n;
        /** The number of interior gridpoints along one edge of the large
         * square. */ 
        const int r;
        /** The total number of interior gridpoints in the large square. */ 
        const int rr;
        /** A status array for the gridpoints. */
        char* const c;
        pss_read(int n_,const char* filename);
        ~pss_read();
        void ascii_art();
        void analyze_glue();
    private:
        /** The x positions of the lower left corners of the component squares.
         */
        std::vector<int> xs;
        /** The y positions of the lower left corners of the component squares.
         */
        std::vector<int> ys;
        /** The side lengths of the squares. */ 
        std::vector<int> ss;
};

#endif
