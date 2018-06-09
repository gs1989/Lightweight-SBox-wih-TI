.func TI_S_Shares
TI_S_:
push {lr}
push {r4-r7}
nop
nop
nop
nop
nop
nop
nop
nop
   movs r3,#0
   mov r4,r0
   eors r3,r4
   mov r4,r0
   movs r6,#8
   rors r4,r6
   mov r5,r0
   movs r6,#16
   rors r5,r6
   ands r4,r5
   eors r3,r4
   mov r4,r0
   movs r6,#24
   rors r4,r6
   mov r5,r1
   ands r4,r5
   eors r3,r4
   mov r4,r0
   mov r5,r0
   movs r6,#16
   rors r5,r6
   ands r4,r5
   eors r3,r4
   mov r4,r0
   movs r6,#8
   rors r4,r6
   mov r5,r0
   movs r6,#24
   rors r5,r6
   ands r4,r5
   eors r3,r4
   mov r4,r0
   movs r6,#16
   rors r4,r6
   mov r5,r1
   ands r4,r5
   eors r3,r4
   mov r4,r0
   movs r6,#24
   rors r4,r6
   mov r5,r1
   movs r6,#8
   rors r5,r6
   ands r4,r5
   eors r3,r4
   mov r4,r0
   mov r5,r0
   movs r6,#24
   rors r5,r6
   ands r4,r5
   eors r3,r4
   mov r4,r0
   movs r6,#16
   rors r4,r6
   mov r5,r1
   movs r6,#8
   rors r5,r6
   ands r4,r5
   eors r3,r4
   mov r4,r0
   movs r6,#8
   rors r4,r6
   mov r5,r1
   movs r6,#16
   rors r5,r6
   ands r4,r5
   eors r3,r4
   mov r4,r0
   mov r5,r1
   movs r6,#16
   rors r5,r6
   ands r4,r5
   eors r3,r4
   mov r4,r0
   movs r6,#8
   rors r4,r6
   mov r5,r1
   movs r6,#24
   rors r5,r6
   ands r4,r5
   eors r3,r4
   mov r4,r0
   mov r5,r1
   movs r6,#24
   rors r5,r6
   ands r4,r5
   eors r3,r4
mov r0,r3
nop
nop
nop
nop
nop
nop
nop
nop
pop {r4-r7}
pop {pc}
.endfunc
Pmask=1010101
.func TI_P_shares
TI_P_shares:
push {lr}
push {r4-r7}
   movs r3,#0 
   movs r2,#1 
   lsls r2,#8 
   subs r2,#1 
   lsls r2,#24 
   ands r2,r0 
   lsrs r2,#24 
   movs r4,#1 
   movs r5,r2 
   ands r5,r4 
   muls r5,r1 
   eors r3,r5 
   lsls r1,#1 
   lsrs r2,#1 
   movs r5,r2 
   ands r5,r4 
   muls r5,r1 
   eors r3,r5 
   lsls r1,#1 
   lsrs r2,#1 
   movs r5,r2 
   ands r5,r4 
   muls r5,r1 
   eors r3,r5 
   lsls r1,#1 
   lsrs r2,#1 
   movs r5,r2 
   ands r5,r4 
   muls r5,r1 
   eors r3,r5 
   lsls r1,#1 
   lsrs r2,#1 
   movs r5,r2 
   ands r5,r4 
   muls r5,r1 
   eors r3,r5 
   lsls r1,#1 
   lsrs r2,#1 
   movs r5,r2 
   ands r5,r4 
   muls r5,r1 
   eors r3,r5 
   lsls r1,#1 
   lsrs r2,#1 
   movs r5,r2 
   ands r5,r4 
   muls r5,r1 
   eors r3,r5 
   lsls r1,#1 
   lsrs r2,#1 
   movs r5,r2 
   ands r5,r4 
   muls r5,r1 
   eors r3,r5 
   lsls r1,#1 
   lsrs r2,#1 
   lsls r0,#8 
   eors r0,r3 
pop {r4-r7}
pop {pc}
.endfunc
