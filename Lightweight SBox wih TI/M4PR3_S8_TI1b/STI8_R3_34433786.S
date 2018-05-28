.func TI_S
TI_S:
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
   str r0,[r1,#0]
   mov r7,r0
   ands r7,r2
   lsrs r7,#1
   mov r4,r0
   ands r4,r3
   lsrs r4,#1
   lsls r5,r4,#3
   eors r5,r4
   eors r7,r5
   nop
   str r7,[r1,#4]
   mov r4,r0
   mov r5,r7
   movs r6,#4
   rors r0,r5
   nop
   str r0,[r1,#8]
   rors r7,r5
   nop
   str r7,[r1,#12]
   rors r0,r5
   nop
   str r0,[r1,#16]
   rors r7,r5
   nop
   str r7,[r1,#20]
   rors r0,r5
   nop
   str r0,[r1,#24]
   rors r7,r5
   nop
   str r7,[r1,#28]
   rors r0,r5
   nop
   str r0,[r1,#32]
   rors r7,r5
   nop
   str r7,[r1,#36]
   rors r0,r5
   nop
   str r0,[r1,#40]
   rors r7,r5
   nop
   str r7,[r1,#44]
   rors r0,r5
   nop
   str r0,[r1,#48]
   rors r7,r5
   nop
   str r7,[r1,#52]
   rors r0,r5
   nop
   str r0,[r1,#56]
   rors r7,r5
   nop
   str r7,[r1,#60]
   movs r0,#0
   nop
   ldr r2,[r1,#8]
   eors r0,r2
   nop
   ldr r2,[r1,#24]
   eors r0,r2
   nop
   ldr r2,[r1,#48]
   eors r0,r2
   nop
   ldr r2,[r1,#4]
   nop
   ldr r3,[r1,#8]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#12]
   nop
   ldr r3,[r1,#16]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#20]
   nop
   ldr r3,[r1,#24]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#44]
   nop
   ldr r3,[r1,#48]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#0]
   nop
   ldr r3,[r1,#8]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#8]
   nop
   ldr r3,[r1,#16]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#16]
   nop
   ldr r3,[r1,#24]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#40]
   nop
   ldr r3,[r1,#48]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#0]
   nop
   ldr r3,[r1,#12]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#4]
   nop
   ldr r3,[r1,#16]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#8]
   nop
   ldr r3,[r1,#20]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#12]
   nop
   ldr r3,[r1,#24]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#16]
   nop
   ldr r3,[r1,#28]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#40]
   nop
   ldr r3,[r1,#52]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#0]
   nop
   ldr r3,[r1,#16]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#8]
   nop
   ldr r3,[r1,#24]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#0]
   nop
   ldr r3,[r1,#20]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#8]
   nop
   ldr r3,[r1,#28]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#12]
   nop
   ldr r3,[r1,#32]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#8]
   nop
   ldr r3,[r1,#32]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#4]
   nop
   ldr r3,[r1,#32]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#8]
   nop
   ldr r3,[r1,#36]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#12]
   nop
   ldr r3,[r1,#40]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#0]
   nop
   ldr r3,[r1,#32]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#8]
   nop
   ldr r3,[r1,#40]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#0]
   nop
   ldr r3,[r1,#36]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#8]
   nop
   ldr r3,[r1,#44]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#20]
   nop
   ldr r3,[r1,#56]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#16]
   nop
   ldr r3,[r1,#56]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#4]
   nop
   ldr r3,[r1,#48]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#12]
   nop
   ldr r3,[r1,#56]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#16]
   nop
   ldr r3,[r1,#60]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#0]
   nop
   ldr r3,[r1,#48]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#8]
   nop
   ldr r3,[r1,#56]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#0]
   nop
   ldr r3,[r1,#52]
   ands r2,r3
   eors r0,r2
   nop
   ldr r2,[r1,#8]
   nop
   ldr r3,[r1,#60]
   ands r2,r3
   eors r0,r2
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
Pmask=88880808
.func TI_P
TI_P:
push {lr}
   movs r5,#0 
   mov r4,r3 
   mov r6,r2 
   ands r4,r0 
   ands r6,r4 
   lsrs r6,#31 
   muls r6,r1 
   lsls r4,#4 
   eors r4,r6 
   eors r5,r4 
   lsrs r1,#1 
   lsrs r2,#1 
   lsrs r3,#1 
   mov r4,r3 
   mov r6,r2 
   ands r4,r0 
   ands r6,r4 
   lsrs r6,#30 
   muls r6,r1 
   lsls r4,#4 
   eors r4,r6 
   eors r5,r4 
   lsrs r1,#1 
   lsrs r2,#1 
   lsrs r3,#1 
   mov r4,r3 
   mov r6,r2 
   ands r4,r0 
   ands r6,r4 
   lsrs r6,#29 
   muls r6,r1 
   lsls r4,#4 
   eors r4,r6 
   eors r5,r4 
   lsrs r1,#1 
   lsrs r2,#1 
   lsrs r3,#1 
   mov r4,r3 
   mov r6,r2 
   ands r4,r0 
   ands r6,r4 
   lsrs r6,#28 
   muls r6,r1 
   lsls r4,#4 
   eors r4,r6 
   eors r5,r4 
   mov r0,r5 
pop {pc}
.endfunc
