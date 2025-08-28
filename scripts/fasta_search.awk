#!/usr/bin/awk -f
BEGIN {
    #print "Searching with regex:", search
    hit = 0
}
{
    if ($0 ~ /^>/) {
        if ($0 ~ search)
        {
            hit = 1
            print $0
        }
        else
	{
            hit = 0
        }
     }
    else if (hit == 1)
    {
        print $0
    }
}

