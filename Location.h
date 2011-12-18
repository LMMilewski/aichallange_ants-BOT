#ifndef LOCATION_H_
#define LOCATION_H_

/*
    struct for representing locations in the grid.
*/
struct Location
{
    int row, col;

    Location()
    {
        row = col = 0;
    };

    Location(int r, int c)
    {
        row = r;
        col = c;
    };

    bool operator< (const Location& loc) const {
        return row < loc.row
            || (row == loc.row && col < loc.col);
    }

    bool operator== (const Location& loc) const {
        return row == loc.row && col == loc.col;
    }

};

inline
std::ostream& operator<< (std::ostream& os, const Location& loc) {
    os << "(" << loc.row << ", " << loc.col << ")";
    return os;
}

// std::ostream& operator<< (std::ostream& os, const Location& loc);

#endif //LOCATION_H_
