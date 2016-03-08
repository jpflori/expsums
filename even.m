function FFInv(t)
    if t eq 0 then
        return 0;
    end if;
    return 1/t;
end function;

function KloostermanSum(a)
    K := Parent(a);
    m := Degree(K);
    E := EllipticCurve([1,0,0,0,a]);
    return #E - 2^m;
end function;

function KloostermanSumPartial(a, b)
    ZZ := Integers();
    K<z> := Parent(a);
    m := Degree(K);
    pow := ZZ!((2^m-1)/3);
    S := 0;
    t := K!1;
    for i in [0..2^m-2] do
        t *:= z;
        if t^pow eq b then
            S +:= (-1)^(ZZ!(Trace(a*t+FFInv(t))));
        end if;
    end for;
    return S;
end function;

function CubicSum(a)
    ZZ := Integers();
    K<z> := Parent(a);
    m := Degree(K);
    S := 0;
    t := K!1;
    for i in [0..2^m-2] do
        t *:= z;
        S +:= (-1)^(ZZ!(Trace(a*t^3)));
    end for;
    return S;
end function;

function CubicSumInv(a, b)
    ZZ := Integers();
    K<z> := Parent(a);
    m := Degree(K);
    S := 0;
    t := K!1;
    for i in [0..2^m-2] do
        t *:= z;
        S +:= (-1)^(ZZ!(Trace(a*t^3 + b/t^3)));
    end for;
    return S;
end function;

//this is the boring sum
//and the best implem?
function S12U1left(a, b, w1, w2)
    ZZ := Integers();
    K := Parent(w1);
    n := Degree(K);
    m1 := ZZ!(n/2);
    m2 := ZZ!(m1/2);
    z := PrimitiveElement(K);
    //pow := ZZ!((2^n-1)/3);
    pow3 := ZZ!((2^m2+1)/3);
    pow2 := 2^m2;
    pow1 := 2^m1-1;
    t1 := z^pow1;
    t1inv := t1^-1;
    t1pow := t1^-2;
    t1powinv := t1^2;
    u1pow := a*w1^2;
    u1powinv := u1pow;
    u1 := K!w2;
    u1inv := K!w2;
    S := 0;
    //sum for u2 \neq 1
    for i in [0..2^(m1-1)-1] do
        u1pow *:= t1pow;
        u1powinv *:= t1powinv;
        u1 *:= t1;
        u1inv *:= t1inv;
        s2 := u1 + u1inv;
//        s2 := b*(s2^-pow);
        s2 := s2^2;
        s2 := s2^pow3;
        s2pow := s2^pow2;
        if s2pow eq b*s2 then
//        if s2 eq K!1 then
            s1 := (-1)^(ZZ!(Trace(u1pow))) + (-1)^(ZZ!(Trace(u1powinv)));
            S +:= s1;
        end if;
    end for;
    return S;
end function;

//this is the boring sum for general extensions
function S12U1leftgen(a, b, K0, Kv)
    ZZ := Integers();
    n := Degree(K0);
    m1 := ZZ!(n/2);
    z := Generator(K0);
    pow := ZZ!((2^n-1)/3);
    t1 := z^(2^m1-1);
    t1pow := t1^(2^m1-1);
    u1 := K0!1;
    u2 := K0!1;
    S := 0;
    //sum for u2 \neq 1
    for i in [0..2^m1-1] do
        u1 *:= t1pow;
        u2 *:= t1;
        s2 := b*Trace(u2,Kv)^(-pow);
        if s2 eq K0!1 then
            s1 := (-1)^(ZZ!(Trace(a*u1)));
            S +:= s1;
        end if;
    end for;
    return S;
end function;

function S12U1partial(a, b, w1, w2);
    ZZ := Integers();
    K<z> := Parent(w1);
    n := Degree(K);
    pow := ZZ!((2^n-1)/3);
    m1 := ZZ!(n/2);
    t1 := z^(2^m1-1);
    t1inv := t1^(-1);
    L<y>, to_K := sub<K|m1>;
    yinv := y^-1;
    w1inv := w1^(-1);
    w2pow := w2^(-pow);
    // sum for u1 \neq w1^-1
    u1 := w1inv;
    u1inv := w1;
    u2 := K!1;
    u2inv := K!1;
    S := 0;
    for i in [0..2^m1-1] do
        u1 *:= t1;
        u1inv *:= t1inv;
        u2 *:= t1;
        u2inv *:= t1inv;
        s2 := b*(u2 + u2inv)^(-pow)*w2pow;
        if s2 eq K!1 then
            s1 := (-1)^(ZZ!(Trace(L!(a*(u1^2+u1inv^2)) + FFInv((L!(u1^2+u1inv^2))))));
            S +:= s1;
        end if;
    end for;
    // sum over F_{2^m}^* \ {1}
    l1 := L!1;
    l1inv := L!1;
    l2 := w1;
    l2inv := w1inv;
    for i in [1..2^m1-2] do
        l1 *:= y;
        l1inv *:= yinv;
        l2 *:= y;
        l2inv *:= yinv;
        s2 := b*(l2 + l2inv)^(-pow)*w2pow;
        if s2 eq K!1 then
            s1 := (-1)^(ZZ!(Trace(L!(a*(l1^2+l1inv^2)) + FFInv((L!(l1^2+l1inv^2))))));
            S +:= s1;
        end if;
    end for;
    return S;
end function;

//this is the boring sum
function SU12(a, b, w1, w2)
    ZZ := Integers();
    K := Parent(w1);
    m0 := Degree(K);
    m1 := ZZ!(m0/2);
    m2 := ZZ!(m1/2);
    L := GF(2, m2);
    M := GF(4);
    z := Generator(K);
    pow := ZZ!((2^m0-1)/3);
    t1 := z^(2^m1-1);
    t1inv := t1^(-1);
    t1pow := t1^(-2);
    t2 := z^((2^m1+1)*(2^m2-1));
    u1 := w1^-1;
    u1pow := u1^(-2);
    u2 := K!1;
    u2inv := K!1;
    S := 0;
    //sum for u2 \neq 1
    for i in [0..2^m1-1] do
        u1pow *:= t1pow;
        u1 *:= t1;
        u2 := K!1;
        for j in [0..2^m2] do
            u2 *:= t2;
            if Trace(u1*u2*w1*w2, L) eq 0 then
                S +:= (-1)^(ZZ!(Trace(a*u1pow))) * (-1)^(ZZ!(Trace(b*M!(u2^pow))));
            end if;
        end for;
    end for;
    return S;
end function;

//this is the boring sum
function SU12trace(a, b, w1, w2)
    ZZ := Integers();
    K := Parent(w1);
    m0 := Degree(K);
    m1 := ZZ!(m0/2);
    m2 := ZZ!(m1/2);
    J := GF(2, m1);
    L := GF(2, m2);
    M := GF(4);
    z := Generator(K);
    pow := ZZ!((2^m0-1)/3);
    t1 := z^(2^m1-1);
    t1inv := t1^(-1);
    t1pow := t1^(-2);
    t2 := z^((2^m1+1)*(2^m2-1));
    u1 := w1^-1;
    u1pow := u1^(-2)*a;
    u1 *:= w1;
    u2 := K!1;
    u2inv := K!1;
    S := 0;
    //sum for u2 \neq 1
    for i in [0..2^m1-1] do
        u1pow *:= t1pow;
        u1 *:= t1;
        S +:= (-1)^(ZZ!(Trace(u1pow))) * (-1)^(ZZ!(Trace(b*M!((w2*Trace(u1,J))^-pow))));
    end for;
    return S;
end function;

// technical stuff to go through cyclotomic classes mod 2^m1
function iatoi(W, n)
    a := 0;
    for i in [1..n] do
        a := ShiftLeft(a,1);
        a +:= W[i];
    end for;
    return a;
end function;

function InitNecklaces(n)
    return [0 : i in [0..n]];
end function;

procedure NextNecklaces(~W, ~a, ~p, n)
    while true do
        p := n;
        while (W[p] eq 1) do
            p := p - 1;
        end while;
        if (p eq 1) then
            W[1] := 1;
            a := iatoi(W, n);
            return;
        end if;
        W[p] := 1;
        for i in [1..(n - p)] do
            W[i + p] := W[i];
        end for;
        if ((n mod p) eq 0) then
            a := iatoi(W, n);
            return;
        end if;
    end while;
end procedure;

// only made for v == 2
function CheckConjecture(v, mv)
    m2 := 2^(v-2)*mv;
    m1 := 2^(v-1)*mv;
    m0 := 2^v*mv;
    ZZ := Integers();
    K1<z1> := GF(2,m1);
    K0<z0> := GF(2,m0);
    w1 := z0^(2^m1-1);
    w2 := z1^(2^m2-1);
    b := GF(4)!1;
    d := -2^m2;
    pow2m := ShiftLeft(1, m1);
    pow2m1 := pow2m - 1;
    W := InitNecklaces(m1);
    j := 0;
    p := 0;
    //for j in [1..2^m1-1] do
    repeat
        NextNecklaces(~W, ~j, ~p, m1);
        if (j mod 3) ne 0 then
            continue;
        end if;
        if (j ne 513) then
            continue;
        end if;
        a := z1^j;
        KS := KloostermanSum(a);
        KSm := -(2^(m2+1)-KS)/3;
        q := (2^(m2-1) - 1)/3;
        if ((KS - 1) mod 3) eq 0 then
            print j, KS;
            w := Cputime();
            for i in [1..2^(m1-1)] do
                Sab := S12U1left(a, b, w1^i, K0!1);
                D := Sab + KSm + (1 - (-1)^(ZZ!Trace(a*w1^(2*i)))) * q;
                if D ne d and D ne 0 then
                    print mv, j, i, D;
                    return false;
                end if;
            end for;
            print Cputime(w);
        end if;
    //end for;
    until (j eq pow2m1);
    return true;
end function;

v := 2;
mv := 7;
m2 := 2^(v-2)*mv;
m1 := 2^(v-1)*mv;
m0 := 2^v*mv;
ZZ := Integers();
Kvm<zvm> := GF(2,mv*2);
K1<z1> := GF(2,m1);
K0<z0> := GF(2,m0);
w1 := z0^(2^m1-1);
w2 := z1^(2^m2-1);

//b := Random(GF(4));
b := GF(4)!1;

/*
a := Random(K1);
KS := KloostermanSum(a);
while ((KS - 1) mod 3) ne 0 do
    a := Random(K1);
    KS := KloostermanSum(a);
end while;
*/

/*
// Looking for the mysterious function
d := -2^m2;
//for j in [0] do
for j in [1..2^m1-1] do
    a := z1^j;
    if (j mod 3) ne 0 then
        continue;
    end if;
    KS := KloostermanSum(a);
    KS0 := KloostermanSum(K0!a);
    KSz1 := KloostermanSum(z1^(ZZ!(j/3)));
    //if KS eq -8 then
    if ((KS - 1) mod 3) eq 0 then
        Sm := {};
        Sp := {};
        print j, a, KS, KS0, KSz1;
        for i in [1..2^(m1-1)] do
            tst := KloostermanSum(w1^i);
            Sab := S12U1left(a, b, w1^(i), K0!1);
            D := Sab - (2^(m2+1)-KS)/3 + (1 - (-1)^(ZZ!Trace(a*w1^(2*i)))) * (2^(m2-1) - 1)/3;
            print j, a, i, KS, KS0, D, tst;//, Trace(z1^(ZZ!(j/3))*w1^(i)), (w1^i + w1^-i)^(ZZ!((2^m1-1)/3)), KloostermanSum(K1!(1/(w1^i+w1^-i)));
            if D eq d then
                Include(~Sm, tst);
            else
                Include(~Sp, tst);
            end if;
        end for;
        print #Sp, #Sm, #(Sp meet Sm);
        print Sp meet Sm;
    end if;
end for;
*/

v := 2;
mv := 3;
m2 := 2^(v-2)*mv;
m1 := 2^(v-1)*mv;
m0 := 2^v*mv;
ZZ := Integers();
Kvm<zvm> := GF(2,mv*2);
K1<z1> := GF(2,m1);
K0<z0> := GF(2,m0);
w1 := z0^(2^m1-1);
w2 := z1^(2^m2-1);
//e := ZZ!((2^m0-1)/3)*(2^m1-1);
e := 2^m1 - 1;
f := ZZ!((2^m0-1)/3);
g := ZZ!((2^m1-1)/3);
for k in [0..2^m1-1] do
    if (k mod 3) ne 0 then
        continue;
    end if;
    a := z1^k;
    c := z1^(ZZ!(k/3));
    KS := KloostermanSum(a);
    if ((KS - 1) mod 3) eq 0 then
        print "===", k/3, a, KS, Trace(1/a), "=================================================";
        for j in [1,2,4,8,16] do//1..2^(m1-1)] do
            print "---", j, Trace(a*w1^(2*j)), Trace(c,GF(4))*GF(4)!(Trace(c*w1^(2*j),K1)^g), "---------";
            for i in [0..2] do //..2] do //2^(m2)] do
                Sab := S12U1left(a, b, w1^(j), w2^(i));
                D := Sab - (2^(m2+1)-KS)/3 + (1 - (-1)^(ZZ!Trace(a*w1^(2*j)))) * (2^(m2-1) - 1)/3;
                print i, D, Trace(a*w1^(2*j)) + Trace(GF(4)!(b*w2^(i*4*ZZ!((2^m2+1)/3))));
            end for;
        end for;
    end if;
end for;
