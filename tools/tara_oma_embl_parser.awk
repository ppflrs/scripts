#Adapted from https://www.biostars.org/p/8241/#8242
/^DE   TARA/ {printf(">%s",$2); next;}
/^(PT|PA)/  {printf(" %s;",$0); next;}
/^\/\// {printf("\n"); next;}
/^    / {printf("\n%s",substr($0,5)); next;}
    {
    /* ignore default */
    }
END   {
    printf("\n");
    }
